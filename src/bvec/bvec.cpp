#include "bvec.h"

// constructor - given a sorted vector of distinct 32bit integers
bvec::bvec(vector<word_t>& vals) {
    count = vals.size();
    // if the density is too low, run length encoding will take MORE space
    if (low_density(vals)) {
        if (DEBUG) printf("constructed non rle\n");
        words = vals;
        rle = false;
        size = words.back();
    }
    else {
        if (DEBUG) printf("construct_rle\n");
        construct_rle(vals);
    }
}

vector<word_t>& bvec::get_words() {
    return words;
}

word_t bvec::get_size() { return size; }
word_t bvec::bytes() { return 4 * words.size(); }

void bvec::compress() {
    if (rle) { /* Throw exception? */ return; }
    vector<word_t> tmp;
    tmp.swap(words);
    construct_rle(tmp);
}

// in place version of the bitwise OR operator.
void bvec::operator|=(bvec& bv) {
    // decide which version we'll be using
    if (rle)
        if (bv.rle)
            rle_OR_rle(bv);
        else
            rle_OR_non(bv);
    else
        if (bv.rle)
            non_OR_rle(bv);
        else
            non_OR_non(bv);
}
// in place version of the bitwise AND operator.
void bvec::operator&=(bvec& bv) {
    // decide which version we'll be using
    if (rle)
        if (bv.rle)
            rle_AND_rle(bv);
        else
            rle_AND_non(bv);
    else
        if (bv.rle)
            non_AND_rle(bv);
        else
            non_AND_non(bv);
}
bvec* bvec::operator|(bvec& rhs) {
    bvec *res = new bvec();
    if (words.size() > rhs.words.size()) {
        res->copy(rhs);
        *res |= *this;
        return res;
    }
    res->copy(*this);
    *res |= rhs;
    return res;
}
bvec* bvec::operator&(bvec& rhs) {
    bvec *res = new bvec();
    if (words.size() > rhs.words.size()) {
        res->copy(rhs);
        *res &= *this;
        return res;
    }
    res->copy(*this);
    *res &= rhs;
    return res;
}

bool bvec::operator==(bvec& other) const {
    return (words == other.words) &&
           (count == other.count) &&
           (size  == other.size)  &&
           (rle   == other.rle);
}

bool bvec::equals(const bvec& other) const {
    return (words == other.words) &&
           (count == other.count) &&
           (size  == other.size)  &&
           (rle   == other.rle);
}

bvec* bvec::copyflip() {
    bvec *res = new bvec();
    res->copy(*this);
    res->flip();
    return res;
}

void bvec::non_OR_non(bvec& bv) {
    
    vector<word_t> res;
    vector<word_t>::iterator a = words.begin();
    vector<word_t>::iterator b = bv.words.begin();
    res.push_back(*a < *b ? *a : *b);
    while(a != words.end() && b != bv.words.end()) {
        if (*a < *b) {
            if (*a != res.back())
                res.push_back(*a);
            ++a;
        }
        else if (*b < *a) {
            if (*b != res.back())
                res.push_back(*b);
            ++b;
        }
        else {
            if (*a != res.back())
                res.push_back(*a);
            ++a;
            ++b;
        }
    }
    if (a != words.end())
        res.insert(res.end(),a,words.end());
    else if (b != bv.words.end())
        res.insert(res.end(),b,bv.words.end());

    count = res.size();
    // TODO: check if it's worth compressing

    words.swap(res);
}

void bvec::non_AND_non(bvec& bv) {
    
    vector<word_t> res;
    vector<word_t>::iterator a = words.begin();
    vector<word_t>::iterator b = bv.words.begin();
    res.push_back(*a < *b ? *a : *b);
    while(a != words.end() && b != bv.words.end()) {
        if (*a < *b)
            ++a;
        else if (*b < *a)
            ++b;
        else {
            res.push_back(*a);
            ++a;
            ++b;
        }
    }
    count = res.size();
    words.swap(res);
}

void bvec::non_AND_rle(bvec& bv) {
    // decompress
    // run non_AND_non
    bvec *tmp = new bvec();
    tmp->copy(bv);
    tmp->decompress();
    non_AND_non(*tmp);
    delete tmp;
}

void bvec::rle_AND_non(bvec& bv) {
    // decompress
    // run non_AND_non
    decompress();
    non_AND_non(bv);
}

void bvec::non_OR_rle(bvec& bv) {
    // compress
    // run rle_OR_rle
    compress();
    rle_OR_rle(bv);
}

void bvec::rle_OR_non(bvec& bv) {
    // compress
    // run rle_OR_rle
    bvec *tmp = new bvec();
    tmp->copy(bv);
    tmp->compress();
    rle_OR_rle(*tmp);
    delete tmp;
}

void bvec::setBit(word_t x) {
    if (rle) { // this should really be done in an rle specific way.
        decompress();
        setBit(x);
        compress();
    }
    else {
        if (words.size() == 0 || x > words.back()) {
            words.push_back(x);
        }
        else {
            vector<word_t>::iterator lb = lower_bound(words.begin(),words.end(),x);
            if (*lb == x) return;
            words.insert(lb,x);
        }
        size++;
        count++;
    }
}

bvec&
bvec::copy(const bvec& bv) {
    words = bv.words;
    count = bv.count;
    size = bv.size;
    rle = bv.rle;
    return *this;
}

ostream& operator<<(ostream &os, const vector<word_t> vec) {
    return os << vec;
}

ostream& operator<<(ostream &os, const bvec &vec) {
    return os << vec;
}

void
bvec::save_to_file(const bvec &bv, const char *filename) {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    oa << bv;
}

void
bvec::restore_from_file(bvec &bv, const char *filename) {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    ia >> bv;
}
