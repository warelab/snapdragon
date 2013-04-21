#include "bvec32.h"
#include "bvec.h"

// constructor - given a previously dumped bvec
bvec::bvec(word_t *buf) {
    word_t nwords = buf[0];
    size = buf[1];
    count = buf[2];
    rle = false;
    if (nwords & BIT1) {
        nwords -= BIT1;
        rle=true;
    }
    words.resize(nwords);
    memcpy(words.data(),buf+3,nwords*4);
    frontier.active_word = words.begin();
    frontier.bit_pos=0;
}

// DIY serialization
size_t
bvec::dump(word_t **buf) {
    // allocate space in buf
    size_t dbytes = sizeof(word_t)*(3 + words.size());
    *buf = (word_t*)malloc(dbytes);
    if (*buf == NULL) {
        fprintf(stderr,"failed to allocate %zi bytes\n",dbytes);
        return 0;
    }
    (*buf)[0] = words.size();
    if (rle) (*buf)[0] |= BIT1;
    (*buf)[1] = size;
    (*buf)[2] = count;
    memcpy(*buf + 3, words.data(), bytes());
    return dbytes;
}

bool
bvec::low_density(vector<word_t>& vals) {
    if (DEBUG) printf("low_density() %u/%u %c %f\n", (word_t)vals.size(),
        vals.back() - vals.front() + 1,
            ((double)vals.size()/(double)(vals.back() - vals.front() + 1) < 1.0/(double)LITERAL_SIZE) ? '<' : '>',
                1.0/(double)LITERAL_SIZE);
    return (double)vals.size()/(double)(vals.back() - vals.front() + 1) < 1.0/(double)LITERAL_SIZE;
    
}

void
bvec::print() {
    printf("rle: %c\n",rle ? 'T' : 'F');
    printf("words:\n");
    for(int i=0;i<words.size();i++) {
        printf(" %i ",i);
        if (rle)
            if ((words[i] & ONEFILL) == ONEFILL) 
                printf("1-fill %zi %zi\n",words[i] & FILLMASK, LITERAL_SIZE*(words[i] & FILLMASK));
            else if (words[i] & BIT1)
                printf("0-fill %zi %zi\n",words[i] & FILLMASK, LITERAL_SIZE*(words[i] & FILLMASK));
            else
                printf("literal %u\n",words[i]);
        else
            printf("direct %u\n",words[i]);
    }
}

void
bvec::construct_rle(vector<word_t>& vals) {
    rle = true;
    frontier.bit_pos=0;
    if (vals.size() == 0) {
        size=0;
        count=0;
        words.clear();
        words.push_back(0); // current word is an empty literal
        frontier.active_word = words.begin();
        return;
    }
    frontier.active_word = words.begin();
    word_t word_end = LITERAL_SIZE - 1;
    word_t word=0;
    word_t gap_words = vals.front()/LITERAL_SIZE;
    if (gap_words > 0) {
        word_end += LITERAL_SIZE*gap_words;
        while (gap_words > FILLMASK) {
            words.push_back(ZEROFULL);
            gap_words -= FILLMASK;
        }
        if (gap_words > 0)
            words.push_back(gap_words | BIT1);
    }
    for(vector<word_t>::iterator ii = vals.begin(); ii != vals.end(); ++ii) {
        if (*ii == word_end)
            word |= 1;
        else if (*ii < word_end)
            word |= ((word_t)1 << (word_end - *ii));
        else {
            if (word == ALL1S)
                if ((words.size() != 0) &&
                    ((words.back() & ONEFILL) == ONEFILL) &&
                    (words.back() != ONEFULL))
                    words.back()++;
                else
                    words.push_back(ONEFILL1);
            else
                words.push_back(word);
            gap_words = (*ii - word_end - 1)/LITERAL_SIZE;
            while (gap_words > FILLMASK) {
                words.push_back(ZEROFULL);
                gap_words -= FILLMASK;
            }
            if (gap_words > 0)
                words.push_back(gap_words | BIT1);
            word_end += (gap_words+1)*LITERAL_SIZE;
            word = (word_end - *ii == LITERAL_SIZE)
                ? 1 : (word_t)1 << (word_end - *ii);
        }
    }
    // add the last word
    if (word == ALL1S) {
        if ((words.size() != 0) &&
            ((words.back() & ONEFILL) == ONEFILL) &&
                (words.back() != ONEFULL))
            words.back()++;
        else
            words.push_back(ONEFILL1);
    } else {
        words.push_back(word);
    }
    size = word_end+1;
}

void
bvec::decompress() {
    if (!rle) { /* Throw exception? */ return; }
    // retrieve the set bits from the compressed vector
    vector<word_t> res;
    res.reserve(cnt());
    word_t pos=0;
    for(vector<word_t>::iterator ii = words.begin(); ii != words.end(); ++ii) {
        if ((*ii & ONEFILL) == ONEFILL) {
            word_t n_ones = LITERAL_SIZE*(*ii & FILLMASK);
            for(word_t i=0;i < n_ones; i++) {
                res.push_back(pos);
                pos++;
            }
        }
        else if (*ii & BIT1)
            pos += LITERAL_SIZE*(*ii & FILLMASK);
        else { // desconstruct the literal word
            for(word_t bit=1; bit<=LITERAL_SIZE; bit++)
                if (*ii & ((word_t)1 << (LITERAL_SIZE-bit)))
                    res.push_back(pos+bit-1);
            pos += LITERAL_SIZE;
        }
    }
    words.swap(res);
    count = words.size();
    rle = false;
}

// only works with compressed - return the position of the next set bit after position x
word_t
bvec::next_one(word_t x) {
    if (!rle) { fprintf(stderr,"next_one() only works on compressed bitvectors\n"); exit(1); }

    x++;
    if (frontier.bit_pos > x) {
        frontier.active_word = words.begin();
        frontier.bit_pos = 0;
    }
    
    while(frontier.active_word != words.end()) {
        // what type of word is it?
        if (*(frontier.active_word) & BIT1) { // fill word
            word_t span = (*(frontier.active_word) & FILLMASK) * LITERAL_SIZE;
            if ((frontier.bit_pos + span >= x) && (*(frontier.active_word) & BIT2)) // x is within a 1-fill
                return x;
            frontier.bit_pos += span;
            x += span;
        }
        else { // literal word
            if ((frontier.bit_pos + LITERAL_SIZE >= x) &&
                ((1UL << (frontier.bit_pos + LITERAL_SIZE - x)) & *(frontier.active_word)))
                return x;
            frontier.bit_pos += LITERAL_SIZE;
            x += LITERAL_SIZE;
        }
        frontier.active_word++;
    }
    return size; // no next set bit, so return the number of bits
}

void
bvec::flip() {
    if (!rle)
        compress();
    for(vector<word_t>::iterator it = words.begin(); it!=words.end(); ++it) {
        if (*it & BIT1) { // fill word - flip BIT2
            if (*it & BIT2) { // 1-fill
                // convert to 0-fill
                *it &= ~BIT2;
            }
            else {
                // convert to 1-fill
                *it |= BIT2;
            }
        }
        else { // literal word - flip all bits
            *it = ~*it;
        }
    }
    count = size-count;
}

void
bvec::matchSize(bvec &bv) {
    if (size < bv.size) {
        word_t gap_words = (bv.size - size)/LITERAL_SIZE;
        while (gap_words > FILLMASK) {
            words.push_back(ZEROFULL);
            gap_words -= FILLMASK;
        }
        if (gap_words > 0)
            words.push_back(BIT1 | gap_words);
        size = bv.size;
    }
    else if (size > bv.size) {
        word_t gap_words = (size - bv.size)/LITERAL_SIZE;
        while (gap_words > FILLMASK) {
            bv.words.push_back(ZEROFULL);
            gap_words -= FILLMASK;
        }
        if (gap_words > 0)
            bv.words.push_back(BIT1 | gap_words);
        bv.size = size;
    }
}

void
bvec::rle_OR_rle(bvec& bv) {
    // ensure that both bvecs are the same size
    this->matchSize(bv);
    if (size == 0)
        return;
    
    vector<word_t> res; // fill this then swap with this.words
    vector<word_t>::iterator a = words.begin();
    vector<word_t>::iterator b = bv.words.begin();

    // maintain the end position of the current word
    word_t a_pos = (*a & BIT1) ? (*a & FILLMASK) : 1;
    word_t b_pos = (*b & BIT1) ? (*b & FILLMASK) : 1;
    word_t res_pos=0;
    word_t next_word;
    bool incr_a = false;
    bool incr_b = false;
    word_t last_pos = size/LITERAL_SIZE;
    while(res_pos != last_pos) {
        if (incr_a) {
            while(a_pos <= res_pos) {
                ++a;
                a_pos += (*a & BIT1) ? *a & FILLMASK : 1;
            }
            incr_a = false;
        }
        if (incr_b) {
            while(b_pos <= res_pos) {
                ++b;
                b_pos += (*b & BIT1) ? *b & FILLMASK : 1;
            }
            incr_b = false;
        }
        if (a_pos == b_pos) {
            if ((*a & ONEFILL) == ONEFILL || (*b & ONEFILL) == ONEFILL)
                next_word = ONEFILL | (a_pos - res_pos);
            else
                if (*a & BIT1)
                    if (*b & BIT1) // zero fill
                        next_word = BIT1 | (a_pos - res_pos);
                    else
                        next_word = *b;
                else
                    if (*b & BIT1) // zero fill
                        next_word = *a;
                    else {
                        word_t u = *a | *b;
                        next_word = (u == ALL1S) ? ONEFILL1 : u;
                    }
            incr_a = true;
            incr_b = true;
            res_pos = a_pos;
        }
        else if (a_pos < b_pos) {
            if ((*b & ONEFILL) == ONEFILL) {
                next_word = ONEFILL | (b_pos - res_pos);
                res_pos = b_pos;
                incr_a = true;
                incr_b = true;
            }
            else { // b is a 0-fill or a literal word
                if ((*a & ONEFILL) == ONEFILL)
                    next_word = ONEFILL | (a_pos - res_pos);
                else if (*a & BIT1) // a is 0-fill
                    next_word = BIT1 | (a_pos - res_pos);
                else // literal
                    next_word = *a;
                res_pos = a_pos;
                incr_a = true;
            }
        }
        else { // a_pos > b_pos
            if ((*a & ONEFILL) == ONEFILL) {
                next_word = ONEFILL | (a_pos - res_pos);
                res_pos = a_pos;
                incr_a = true;
                incr_b = true;
            }
            else { // a is a 0-fill or a literal word
                if ((*b & ONEFILL) == ONEFILL)
                    next_word = ONEFILL | (b_pos - res_pos);
                else if (*b & BIT1) // b is 0-fill
                    next_word = BIT1 | (b_pos - res_pos);
                else // literal
                    next_word = *b;
                res_pos = b_pos;
                incr_b = true;
            }
        }
        if ((next_word & BIT1)
            && res.size() > 0
            && (res.back() & BIT1)
            && (res.back() & BIT2) == (next_word & BIT2)
            && (res.back() & FILLMASK) < FILLMASK
        ) {
            // merge fill words
            // but just don't exceed the capacity
            word_t n_words = (res.back() & FILLMASK) + (next_word & FILLMASK);
            if (n_words >= FILLMASK) {
                res.back() |= FILLMASK;
                n_words -= FILLMASK;
                if (n_words > 0)
                     res.push_back((next_word & ONEFILL) | n_words);
            }
            else
                res.push_back(next_word);
        }
        else
            res.push_back(next_word);
    }
    words.swap(res);
    count=0;
    
    // decide whether to decompress
}

// in place version of the bitwise AND operator.
void
bvec::rle_AND_rle(bvec& bv) {
    // ensure that both bvecs are the same size
    this->matchSize(bv);
    if (size == 0)
        return;
    
    vector<word_t> res; // fill this then swap with this.words
    vector<word_t>::iterator a = words.begin();
    vector<word_t>::iterator b = bv.words.begin();

    // maintain the end position of the current word
    word_t a_pos = (*a & BIT1) ? *a & FILLMASK : 1;
    word_t b_pos = (*b & BIT1) ? *b & FILLMASK : 1;
    word_t res_pos=0;
    word_t next_word;
    bool incr_a = false;
    bool incr_b = false;
    word_t last_pos = size/LITERAL_SIZE;
    while(res_pos != last_pos) {
        if (incr_a) {
            while(a_pos <= res_pos) {
                ++a;
                a_pos += (*a & BIT1) ? *a & FILLMASK : 1;
            }
            incr_a = false;
        }
        if (incr_b) {
            while(b_pos <= res_pos) {
                ++b;
                b_pos += (*b & BIT1) ? *b & FILLMASK : 1;
            }
            incr_b = false;
        }
        if (a_pos == b_pos) {
            if ((*a & ONEFILL) == ONEFILL && (*b & ONEFILL) == ONEFILL)
                next_word = ONEFILL | (a_pos - res_pos);
            else if ((*a & BIT1) || (*b & BIT1))
                next_word = BIT1 | (a_pos - res_pos);
            else {
                word_t u = *a & *b;
                next_word = (u == 0) ? BIT1 | 1 : u;
            }
            incr_a = true;
            incr_b = true;
            res_pos = a_pos;
        }
        else if (a_pos < b_pos) {
            if ((*b & ONEFILL) == ONEFILL) {
                if((*a & ONEFILL) == ONEFILL)
                    next_word = ONEFILL | (a_pos - res_pos);
                else if (*a & BIT1)
                    next_word = BIT1 | (a_pos - res_pos);
                else
                    next_word = *a;
                res_pos = a_pos;
                incr_a = true;
            }
            else { // b is a 0-fill because it can't be a literal word and have b_pos > a_pos
                next_word = BIT1 | (b_pos - res_pos);
                res_pos = b_pos;
                incr_a = true;
                incr_b = true;
            }
        }
        else { // a_pos > b_pos
            if ((*a & ONEFILL) == ONEFILL) {
                if((*b & ONEFILL) == ONEFILL)
                    next_word = ONEFILL | (b_pos - res_pos);
                else if (*b & BIT1)
                    next_word = BIT1 | (b_pos - res_pos);
                else
                    next_word = *b;
                res_pos = b_pos;
                incr_b = true;
            }
            else { // a is a 0-fill because it can't be a literal word and have a_pos > b_pos
                next_word = BIT1 | (a_pos - res_pos);
                res_pos = a_pos;
                incr_a = true;
                incr_b = true;
            }
        }
        if ((next_word & BIT1)
            && res.size() > 0
            && (res.back() & BIT1)
            && (res.back() & BIT2) == (next_word & BIT2)
            && (res.back() & FILLMASK) < FILLMASK
        ) {
            // merge fill words
            // but just don't exceed the capacity
            word_t n_words = (res.back() & FILLMASK) + (next_word & FILLMASK);
            if (n_words >= FILLMASK) {
                res.back() |= FILLMASK;
                n_words -= FILLMASK;
                if (n_words > 0)
                     res.push_back((next_word & ONEFILL) | n_words);
            }
            else
                res.push_back(next_word);
        }
        else
            res.push_back(next_word);
    }
    words.swap(res);
    count=0;
}

bool
bvec::find(word_t x) {
    if (!rle) return binary_search(words.begin(), words.end(), x);
    // This function may be called on a sequence of increasing values
    // Use a checkpoint to determine if we can start at the active word
    // or if we need to go back to words.begin().
    if (frontier.bit_pos > x) {
        frontier.active_word = words.begin();
        frontier.bit_pos = 0;
    }
    while(frontier.active_word != words.end()) {
        // what type of word is it?
        if (*(frontier.active_word) & BIT1) { // fill word
            word_t span = (*(frontier.active_word) & FILLMASK) * LITERAL_SIZE;
            if (frontier.bit_pos + span >= x) {
                if (*(frontier.active_word) & BIT2) // 1-fill
                    return true;
                return false;
            }
            frontier.bit_pos += span;
        }
        else { // literal word
            if (frontier.bit_pos + LITERAL_SIZE >= x) {
                if ((1UL << (frontier.bit_pos + LITERAL_SIZE - x)) & *(frontier.active_word))
                    return true;
                return false;
            }
            frontier.bit_pos += LITERAL_SIZE;
        }
        frontier.active_word++;
    }
    return false;
}

// use the frontier instead of words.back() because we have frontier.bit_pos
// advance the frontier after appending
void
bvec::appendFill(bool bit, word_t n) {
    if (bit) count += n;
    if ((*(frontier.active_word) & BIT1) == 0) { // current word is a literal word
        word_t bits_available = LITERAL_SIZE - frontier.bit_pos; // size % LITERAL_SIZE;
//      fprintf(stderr,"bits_available: %u\n",bits_available);
        if(bit) {
            // append some ones
            if (n < bits_available) {
                size += n;
                frontier.bit_pos += n;
                *(frontier.active_word) |= ((1UL << n) - 1) << (bits_available - n);
                return;
            }
            else {
                size += bits_available;
                n -= bits_available;
                frontier.bit_pos = 0;
                *(frontier.active_word) |= (1 << bits_available) - 1;
                // the word is full, check if we should convert to a 1-fill
                if (*(frontier.active_word) == ALL1S) {
                    *(frontier.active_word) = ONEFILL1;
                }
            }
        }
        else {
            // append some zeros
            if (n < bits_available) {
                size += n;
                frontier.bit_pos += n;
                return;
            }
            else {
                size += bits_available;
                n -= bits_available;
                frontier.bit_pos = 0;
                // the word is full, check if we should convert to a 0-fill
                if (*(frontier.active_word) == 0) {
                    *(frontier.active_word) = ZEROFILL1;
                }
            }
        }
    }
//  print();
    if (n==0) return;
    // append/update fill words
    word_t n_fills = n/LITERAL_SIZE;
    if (n_fills > 0) {
        if (bit) {
            if (*(frontier.active_word) & ONEFILL) { // extend previous 1-fill
                *(frontier.active_word) += n_fills;
            }
            else { // append a 1-fill
                words.push_back(ONEFILL | n_fills);
                frontier.active_word = words.end() - 1;
            }
        }
        else {
            if (*(frontier.active_word) & BIT1 && !(*frontier.active_word & BIT2)) { // extend previous 0-fill
                *(frontier.active_word) += n_fills;
            }
            else { // append a 1-fill
                words.push_back(BIT1 | n_fills);
                frontier.active_word = words.end()-1;
            }
        }
        n -= n_fills*LITERAL_SIZE;
        size += n_fills*LITERAL_SIZE;
    }
    // add the remaining bits to a literal word
    if (n>0) {
        if (bit)
            words.push_back(((1<<n)-1) << (LITERAL_SIZE - n));
        else
            words.push_back(0);
        size += n;
        frontier.active_word = words.end() - 1;
        frontier.bit_pos = n;
    }
}

word_t bvec::cnt() {
   if (count == 0)
       for(vector<word_t>::iterator it = words.begin(); it != words.end(); ++it)
           count += (*it & BIT2) ? (*it & BIT1) ? (*it & FILLMASK) * LITERAL_SIZE : 0 : __builtin_popcount(*it);
   return count;
}
