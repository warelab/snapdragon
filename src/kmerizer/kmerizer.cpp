#include "kmerizer.h"
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <bitset>
#include <algorithm> // stl sort
#include <cstdio>  // sprintf()
#include <cstring> // memcpy()
#include <sys/stat.h> // mkdir()
#include <sys/time.h> // gettimeofday()

Kmerizer::Kmerizer(const size_t klength,
                   const size_t threads,
                   const char * outdir,
                   const char   mode)
{
    this->k       = klength-4;
    this->outdir  = new char[strlen(outdir)]; strcpy(this->outdir, outdir);
    this->mode    = mode;
    this->threads = threads;

    this->kmask       = 0xFFFFFFFFFFFFFFFFULL;
    this->shiftlastby = 62;
    this->lshift      = 0;
    this->rshift      = 0;
    if (k % 32 > 0) {
        kmask = (1ULL << (2 * (k % 32))) - 1;
        shiftlastby = 2 * (k % 32) - 2;
        lshift = 2 * (k % 32);
        rshift = 64 - lshift;
    }
    this->nwords     = ((k-1)>>5)+1;
    this->kmerSize   = this->nwords * sizeof(kword_t);
    this->threadBins = NBINS / this->threads;
    this->state      = READ;
    
    fprintf(stderr,"new() k: %zi nwords: %zi, kmerSize: %zi\n",k,nwords,kmerSize);
}

int Kmerizer::allocate(const size_t maximem) {
    fprintf(stderr, "Kmerizer::allocate(%zi)", maximem);
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);

    // allocate memory for the common kmer lookup table.
    // allocate memory for the kmer buffers - assume uniform distribution
    // zero the tallies
    memset(binTally, 0, sizeof(uint32_t) * NBINS);
    memset(rareKmers, 0, sizeof(uint32_t) * NBINS);
    memset(lutTally, 0, sizeof(uint32_t) * NBINS);
    memset(batches, 0, sizeof(size_t) * NBINS);
    uint32_t capacity = maximem / kmerSize / NBINS;
    for (size_t i = 0; i < NBINS; i++) {
        binCapacity[i] = capacity - (capacity >> 2);
        totalCapacity[i] = capacity;
        kmerBuf[i] = (kword_t *) calloc(binCapacity[i], kmerSize);
        if (kmerBuf[i] == NULL) return 1;
        kmerLutK[i] = (kword_t *) calloc(binCapacity[i], kmerSize);
        if (kmerLutK[i] == NULL) return 1;
        kmerLutV[i] = (uint32_t *) calloc(binCapacity[i], sizeof(uint32_t));
        if (kmerLutV[i] == NULL) return 1;
        
        lutTally[i] = 1; // heuristic: add the A* kmer in each bin
    }

    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    fprintf(stderr," took %f seconds\n",elapsedTime);
    return 0;
}

kword_t * Kmerizer::canonicalize(kword_t *packed, kword_t *rcpack, size_t *bin) const {

    // copy the words (reverse their order)
    for (size_t i=0;i<nwords;i++)
        rcpack[i] = packed[nwords-1-i];

    // shift bits from word to word if necessary
    if (lshift) {
        for (size_t i=0;i<nwords-1;i++) {
            rcpack[i] |= rcpack[i+1] << lshift;
            rcpack[i+1] >>= rshift;
        }
    }
    // bitwise reverse complement of values from 0 to 255
    static const kword_t rctable[256] = {
      255,191,127,63,239,175,111,47,223,159,95,31,207,143,79,15,
      251,187,123,59,235,171,107,43,219,155,91,27,203,139,75,11,
      247,183,119,55,231,167,103,39,215,151,87,23,199,135,71,7,
      243,179,115,51,227,163,99,35,211,147,83,19,195,131,67,3,
      254,190,126,62,238,174,110,46,222,158,94,30,206,142,78,14,
      250,186,122,58,234,170,106,42,218,154,90,26,202,138,74,10,
      246,182,118,54,230,166,102,38,214,150,86,22,198,134,70,6,
      242,178,114,50,226,162,98,34,210,146,82,18,194,130,66,2,
      253,189,125,61,237,173,109,45,221,157,93,29,205,141,77,13,
      249,185,121,57,233,169,105,41,217,153,89,25,201,137,73,9,
      245,181,117,53,229,165,101,37,213,149,85,21,197,133,69,5,
      241,177,113,49,225,161,97,33,209,145,81,17,193,129,65,1,
      252,188,124,60,236,172,108,44,220,156,92,28,204,140,76,12,
      248,184,120,56,232,168,104,40,216,152,88,24,200,136,72,8,
      244,180,116,52,228,164,100,36,212,148,84,20,196,132,68,4,
      240,176,112,48,224,160,96,32,208,144,80,16,192,128,64,0,
    };

    size_t rcbin = rctable[rcpack[0] & 255UL];
    if (rcbin > *bin)
        return packed;

    // shift each word by 8 bits
    rcpack[0] >>= 8;
    for (size_t i=0;i<nwords-1;i++) {
        rcpack[i] |= rcpack[i+1] << 56;
        rcpack[i+1] >>= 8;
    }
    if (lshift >= 8)
        rcpack[nwords-1] |= *bin << (lshift-8);
    else if (lshift > 0) {
        rcpack[nwords-2] |= *bin << 56+lshift;
        rcpack[nwords-1] |= *bin >> lshift;
    }
    else
        rcpack[nwords-1] |= *bin << 56;

    int cmp=0;
    if (rcbin > *bin)
        cmp = 1;
    for (size_t i=0;i<nwords;i++) {
        rcpack[i] = revcomp(rcpack[i]);
        if (i==nwords-1) rcpack[i] >>=rshift;
        if (cmp == 0) {
            if (packed[i] < rcpack[i])
                cmp = -1;
            else if (packed[i] > rcpack[i])
                cmp = 1;
        }
    }
    if(cmp > 0) {
        *bin = rcbin;
        return rcpack;
    }
    else
        return packed;
}

size_t Kmerizer::nextKmer(kword_t* kmer, size_t bin, const char nucl) {
    // update the bin
    bin <<= 2;
    if (nwords==1)
        bin |= (kmer[0] >> shiftlastby);
    else
        bin |= (kmer[0] >> 62);
    bin &= 255; // mask out the higher order bits
    
    // shift first kmer by 2 bits
    kmer[0] <<= 2;
    for (size_t w = 1; w < nwords - 1; w++) { // middle (full length) words
        kmer[w-1] |= (kmer[w] >> 62); // move the top two bits
        kmer[w] <<= 2; // make room for the next word
    }
    // last word (if not the first)
    if(nwords>1) {
        kmer[nwords-2] |= (kmer[nwords-1] >> shiftlastby);
        kmer[nwords-1] <<= 2;
    }
    kmer[nwords-1] |= twoBit(nucl);
    kmer[nwords-1] &= kmask;
    return bin;
}

void Kmerizer::addSequence(const char* seq, const int length) {
    if (length < k+4) return;
    kword_t packed[nwords];
    kword_t rcpack[nwords];
    memset(packed, 0, kmerSize);
    size_t bin=0;
    for (size_t i = 0; i < k+4; i++)
        bin = nextKmer(packed,bin,seq[i]);

    kword_t *kmer = packed;
    if (mode == CANONICAL)
        kmer = canonicalize(packed,rcpack,&bin);

    memcpy(kmerBuf[bin] + nwords*binTally[bin], kmer, kmerSize);
    binTally[bin]++;
    if (binTally[bin] == binCapacity[bin]) uniqify();
    // pack the rest of the sequence
    for (size_t i=k+4; i<(size_t)length;i++) {
        bin = nextKmer(packed, bin, seq[i]);
        if (mode == CANONICAL)
            kmer = canonicalize(packed,rcpack, &bin);
        memcpy(kmerBuf[bin] + nwords*binTally[bin], kmer, kmerSize);
        binTally[bin]++;
        if (binTally[bin] == binCapacity[bin]) uniqify();
    }
}

int kmercmp(const void *k1, const void *k2, size_t nwords) {
    for (size_t i=0;i<nwords;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}

int Kmerizer::searchLut(kword_t *kmer, size_t bin) {
    kword_t* A = kmerLutK[bin];
    int imin = 0;
    int imax = lutTally[bin];
    while (imin < imax) {
        int imid = (imin+imax)>>1;
        int cmp = kmercmp(A + nwords*imid, kmer, nwords);
        if (cmp < 0)
            imin = imid + 1;
        else if (cmp == 0)
            return imid;
        else
            imax = imid;
    }
    return -1;
}


void Kmerizer::save() {
    uniqify();
    mergeBatches();
}

void Kmerizer::load() {
    fprintf(stderr,"Kmerizer::load()");
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doLoadIndex, this, i, j));
    }
    tg.join_all();
    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    fprintf(stderr," took %f seconds\n",elapsedTime);
}

// given one kmer, pack it, canonicalize it, hash it, find it
/*
uint32_t Kmerizer::find(const char* seq) {
    // pack it
    kword_t packed[nwords];
    kword_t rcpack[nwords];
    memset(packed,0,kmerSize);
    for (size_t i = 0; i < k; i++)
        nextKmer(packed,seq[i]);
    kword_t *kmer = packed;

    // canonicalize it
    if (mode == CANONICAL) {
        kmer = canonicalize(packed,rcpack);
    }

    // hash it
    size_t bin = hashkmer(kmer,0);

    // check if we've loaded the bit slices for this bin
    if (slices[bin].empty()) {
        // load the kmer slices for this bin
        char fname[100];
        sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
        vector<uint32_t> junk;
        readBitmap(fname,junk,slices[bin]);
    }

    // search the index by iterative boolean operations
    BitVector *res = new BitVector(true); // first create an empty bitvector
    res->appendFill(true,slices[bin][0]->getSize()); // then make it all 1's
    for (size_t w=0;w<nwords;w++) {
        for (size_t b=0;b<64;b++) {
            if (kmer[w] & (1ULL << 63-b))
                *res &= *(slices[bin][w*nwords + b]);
            else
                *res &= *(slices[bin][w*nwords + b]->copyflip());
            if (res->cnt() == 0) return 0;
        }
    }
    // if we've reached this point, there should be one set bit in res
    if (res->cnt() != 1) {
        fprintf(stderr,"expected only 1, but cnt() returned %u\n",res->cnt());
        exit(1);
    }
    // need to figure out which bit
    // and lookup the associated count
    res->decompress();
    vector<uint32_t> hits = res->getWords();
    return frequency(bin,hits[0]);
}
*/

uint32_t Kmerizer::frequency(size_t bin, uint32_t pos) {
    size_t key=1;
    while (key<kmerFreq[bin].size())
        if (counts[bin][key]->find(pos))
            key++;
        else
            break;
    return kmerFreq[bin][key-1];
}


void Kmerizer::histogram() {

    // merge the NBINS counts bvecs
    uint32_t todo = NBINS;
    size_t offset[NBINS];
    size_t done[NBINS];
    memset(offset,0,sizeof(size_t)*NBINS);
    memset(done,0,sizeof(size_t)*NBINS);
    uint32_t key=1; // min kmer frequency
    while (todo > 0) {
        uint32_t val = 0;
        for (size_t i=0;i<NBINS;i++) {
            if (done[i]==0) {
                if (kmerFreq[i][offset[i]] == key) {
                    val += counts[i][offset[i]]->cnt();
                    if (offset[i]<kmerFreq[i].size()-1) // undo the range encoding <=
                        val -= counts[i][offset[i]+1]->cnt();
                    offset[i]++;
                    if (offset[i] == kmerFreq[i].size()) {
                        done[i]=1;
                        todo--;
                    }
                }
            }
        }
        if (val > 0) {
            printf("%u %u\n",key,val);
        }
        key++;
    }
}

void Kmerizer::serialize() {
    // unique the batch
    uniqify();
    // write to disk
    writeBatch();

    // zero the data
    memset(binTally,0,sizeof(uint32_t)*NBINS);
    for (size_t i=0;i<NBINS;i++) {
        for (size_t j=0; j<counts[i].size(); j++)
            delete counts[i][j]; // BitVector destructor
    }
    state = READ;
}

void Kmerizer::uniqify() {
//    fprintf(stderr,"Kmerizer::uniqify()");
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);
    boost::thread_group tg;

    // for (size_t i = 0; i < NBINS; i += threadBins) {
    //     size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
    //     // fprintf(stderr,"i: %zi, threadBins: %zi, j: %zi\n",i,threadBins,j);
    //     tg.create_thread(boost::bind(&Kmerizer::doUnique, this, i, j));
    // }
    for(size_t i=0;i<NBINS;i++)
        if (binTally[i] > 0) //0.5*(float)binCapacity[i])
            tg.create_thread(boost::bind(&Kmerizer::doUnique, this, i,i+1));
    tg.join_all();
    state = QUERY;
    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
 //   fprintf(stderr," took %f seconds\n",elapsedTime);
}

void Kmerizer::writeBatch() {
    fprintf(stderr,"Kmerizer::writeBatch()");
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doWriteBatch, this, i, j));
    }
    tg.join_all();
    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    fprintf(stderr," took %f seconds\n",elapsedTime);
}

void Kmerizer::mergeBatches() {
    fprintf(stderr,"Kmerizer::mergeBatches()");
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);
    boost::thread_group tg;
    for (size_t i=0;i<NBINS;i+=threadBins) {
        size_t j=(i+threadBins>NBINS) ? NBINS : i+threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doMergeBatches, this, i, j));
    }
    tg.join_all();
    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    fprintf(stderr," took %f seconds\n",elapsedTime);
}


int compare_kmers1(const void *k1, const void *k2) {
    if (*(kword_t*)k1 > *(kword_t*)k2) return 1;
    if (*(kword_t*)k1 < *(kword_t*)k2) return -1;
    return 0;
}
int compare_kmers2(const void *k1, const void *k2) {
    for (size_t i=0;i<2;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers3(const void *k1, const void *k2) {
    for (size_t i=0;i<3;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers4(const void *k1, const void *k2) {
    for (size_t i=0;i<4;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers5(const void *k1, const void *k2) {
    for (size_t i=0;i<5;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers6(const void *k1, const void *k2) {
    for (size_t i=0;i<6;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers7(const void *k1, const void *k2) {
    for (size_t i=0;i<7;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}
int compare_kmers8(const void *k1, const void *k2) {
    for (size_t i=0;i<8;i++) {
        if (*((kword_t*)k1+i) < *((kword_t*)k2+i)) return -1;
        if (*((kword_t*)k1+i) > *((kword_t*)k2+i)) return 1;
    }
    return 0;
}


struct kmer1_t { kword_t w[1]; };
struct kmer2_t { kword_t w[2]; };
struct kmer3_t { kword_t w[3]; };
struct kmer4_t { kword_t w[4]; };
struct kmer5_t { kword_t w[5]; };
struct kmer6_t { kword_t w[6]; };
struct kmer7_t { kword_t w[7]; };
struct kmer8_t { kword_t w[8]; };

bool operator<(const kmer1_t& a, const kmer1_t& b) {
    return a.w[0] < b.w[0];
}

bool operator<(const kmer2_t& a, const kmer2_t& b) {
    for (size_t i=0;i<2;i++) {
        if (a.w[i] < b.w[i]) return true;
        if (a.w[i] > b.w[i]) return false;
    }
    return false;
}
bool operator<(const kmer3_t& a, const kmer3_t& b) {
    for (size_t i=0;i<3;i++) {
        if (*(&a.w+i) < *(&b.w+i)) return true;
        if (*(&a.w+i) > *(&b.w+i)) return false;
    }
    return false;
}
bool operator<(const kmer4_t& a, const kmer4_t& b) {
    for (size_t i=0;i<4;i++) {
        if (*(&a.w+i) < *(&b.w+i)) return true;
        if (*(&a.w+i) > *(&b.w+i)) return false;
    }
    return false;
}
bool operator<(const kmer5_t& a, const kmer5_t& b) {
    for (size_t i=0;i<5;i++) {
        if (*(&a.w+i) < *(&b.w+i)) return true;
        if (*(&a.w+i) > *(&b.w+i)) return false;
    }
    return false;
}
bool operator<(const kmer6_t& a, const kmer6_t& b) {
    for (size_t i=0;i<6;i++) {
        if (*(&a.w+i) < *(&b.w+i)) return true;
        if (*(&a.w+i) > *(&b.w+i)) return false;
    }
    return false;
}
bool operator<(const kmer7_t& a, const kmer7_t& b) {
    for (size_t i=0;i<7;i++) {
        if (*(&a.w+i) < *(&b.w+i)) return true;
        if (*(&a.w+i) > *(&b.w+i)) return false;
    }
    return false;
}
bool operator<(const kmer8_t& a, const kmer8_t& b) {
    for (size_t i=0;i<8;i++) {
        if (*(&a.w+i) < *(&b.w+i)) return true;
        if (*(&a.w+i) > *(&b.w+i)) return false;
    }
    return false;
}

void Kmerizer::doUnique(const size_t from, const size_t to) {
//    fprintf(stderr,"doUnique(%zi,%zi)\n",from, to);
    for (size_t bin=from; bin<to; bin++) {
//        fprintf(stderr,"nwords: %zi, binTally[%zi]: %u, ztally: %u\n",nwords,bin,binTally[bin],kmerLutV[bin][0]);
        
        // filter the kmers that are in the LUT
        // shift kmers that are not in the LUT to the beginning of the buffer so we can sort them
        kword_t *kmer;
        if (lutTally[bin]==1) { // first time, check if a kmer is 0
            if (nwords==1) {
                for (size_t i=rareKmers[bin]; i<binTally[bin]; i++) {
                    kmer = kmerBuf[bin] + i;
                    if (*kmer != 0) {
                        memcpy(kmerBuf[bin] + nwords*rareKmers[bin], kmer, kmerSize);
                        rareKmers[bin]++;
                    }
                    else
                        kmerLutV[bin][0]++;
                }
            }
            else {
                for (size_t i=rareKmers[bin]; i<binTally[bin]; i++) {
                    kmer = kmerBuf[bin] + i*nwords;
                    int rc=0;
                    for(int i=0;i<nwords;i++) {
                        if (*(kmer+i) != 0) {
                            rc = -1;
                            break;
                        }
                    }
                    if (rc < 0) {
                        memcpy(kmerBuf[bin] + nwords*rareKmers[bin], kmer, kmerSize);
                        rareKmers[bin]++;
                    }
                    else
                        kmerLutV[bin][rc]++;
                }
            }
        }
        else {
            if (nwords==1) {
                for (size_t i=rareKmers[bin]; i<binTally[bin]; i++) { // if it's not the first time doing this we don't have to check ALL kmers
                    kmer = kmerBuf[bin] + i*nwords;
                    int rc = searchLut1(kmer,bin);
                    if (rc < 0) { // didn't find it
                        // maybe instead of copying them, just keep track of the gaps so they can be filled later
                        memcpy(kmerBuf[bin] + nwords*rareKmers[bin], kmer, kmerSize);
                        rareKmers[bin]++;
                    }
                    else
                        kmerLutV[bin][rc]++;
                }
            }
            else {
                for (size_t i=rareKmers[bin]; i<binTally[bin]; i++) { // if it's not the first time doing this we don't have to check ALL kmers
                    kmer = kmerBuf[bin] + i*nwords;
                    int rc = searchLut(kmer,bin);
                    if (rc < 0) { // didn't find it
                        // maybe instead of copying them, just keep track of the gaps so they can be filled later
                        memcpy(kmerBuf[bin] + nwords*rareKmers[bin], kmer, kmerSize);
                        rareKmers[bin]++;
                    }
                    else
                        kmerLutV[bin][rc]++;
                }
                
            }
        }
        binTally[bin] = rareKmers[bin];
        // need to decide whether to continue or return to filling the buffer
        // compare binTally[bin] to binCapacity[bin] - what ratio? 90%?
        if (binTally[bin] >= 0.9*(float)binCapacity[bin]) {
            rareKmers[bin]=0;
            // if we have gaps in the kmerBuf[bin], we need to memmove the blocks into place before sorting
//            fprintf(stderr,"binTally[%zi]: %u > %f\n",bin,binTally[bin],0.9*(float)binCapacity[bin]);
            if (1) { // stl sort
                if( nwords == 1 ) sort((kmer1_t*)kmerBuf[bin], (kmer1_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 2 ) sort((kmer2_t*)kmerBuf[bin], (kmer2_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 3 ) sort((kmer3_t*)kmerBuf[bin], (kmer3_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 4 ) sort((kmer4_t*)kmerBuf[bin], (kmer4_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 5 ) sort((kmer5_t*)kmerBuf[bin], (kmer5_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 6 ) sort((kmer6_t*)kmerBuf[bin], (kmer6_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 7 ) sort((kmer7_t*)kmerBuf[bin], (kmer7_t*)(kmerBuf[bin] + nwords * binTally[bin]));
                if( nwords == 8 ) sort((kmer8_t*)kmerBuf[bin], (kmer8_t*)(kmerBuf[bin] + nwords * binTally[bin]));
            }
            if (0) { // qsort
                if( nwords == 1 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers1);
                if( nwords == 2 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers2);
                if( nwords == 3 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers3);
                if( nwords == 4 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers4);
                if( nwords == 5 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers5);
                if( nwords == 6 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers6);
                if( nwords == 7 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers7);
                if( nwords == 8 ) qsort(kmerBuf[bin], binTally[bin], kmerSize, compare_kmers8);
            }
            // after the first batch, add the kmers to the common kmer lookup table
            if (lutTally[bin] == 1) {
                // sampling phase is done, release constraint on binCapacity (for when there's lots of memory)
                // uniq -c
                // instead of selecting the most frequent kmers for the lookup table, just take them all
                // the sampling phase will not use excessive RAM, (re)allocate when finished
                uint32_t distinct = 0;
                for (size_t i=1;i<binTally[bin];i++) {
                    kword_t *ith = kmerBuf[bin] + i*nwords;
                    if(kmercmp(kmerLutK[bin] + distinct*nwords, ith, nwords) == 0)
                        kmerLutV[bin][distinct]++;
                    else {
                        distinct++;
                        memcpy(kmerLutK[bin] + distinct*nwords, ith, kmerSize);
                        kmerLutV[bin][distinct] = 1;
                    }
                }
                lutTally[bin] = distinct+1;
                // resize the kmerLutK and kmerLutV arrays to hold lutTally[bin] items.
                kmerLutK[bin] = (kword_t*) realloc(kmerLutK[bin], lutTally[bin]*kmerSize);
                if (kmerLutK[bin] == NULL) {
                    // error in realloc
                    fprintf(stderr,"error in kmerLutK realloc");
                    exit(3);
                }
                kmerLutV[bin] = (uint32_t*) realloc(kmerLutV[bin], lutTally[bin]*sizeof(uint32_t));
                if (kmerLutV[bin] == NULL) {
                    // error in realloc
                    fprintf(stderr,"error in kmerLutV realloc");
                    exit(3);
                }
                binCapacity[bin] = totalCapacity[bin] - lutTally[bin];
                // reallocate kmerBuf[bin] to hold binCapacity[bin] kmers
//                fprintf(stderr,"new binCapacity[%zi]=%u\n",bin,binCapacity[bin]);
                kmerBuf[bin] = (kword_t*) realloc(kmerBuf[bin], binCapacity[bin] * kmerSize);
                if (kmerBuf[bin] == NULL) {
                    // error in realloc
                    fprintf(stderr,"error in kmerBuf realloc");
                    exit(3);
                }
            }
            else { // not in the sampling phase, so kmers are "rare". use 8 bits for the tally
                // uniq -c
                uint32_t distinct = 0;
                vector<uint8_t> tally;
                tally.push_back(1); // first kmer
                for (size_t i=1;i<binTally[bin];i++) {
                    kword_t *ith = kmerBuf[bin] + i*nwords;
                    if(tally.back() < 255 && kmercmp(kmerBuf[bin] + distinct*nwords, ith, nwords) == 0)
                        tally.back()++;
                    else {
                        distinct++;
                        tally.push_back(1);
                        memcpy(kmerBuf[bin] + distinct*nwords, ith, kmerSize);
                    }
                }
                distinct++;
                binTally[bin] = distinct;

                if (binTally[bin] > 0) { // its possible that all values went into the kmerLut
                    batches[bin]++;

                    FILE *fp;

                    // open output file for kmerBuf
                    char kmer_file[100];
                    sprintf(kmer_file,"%s/%zi-mers.%zi.%zi.raw",outdir,k,bin,batches[bin]);
                    fp = fopen(kmer_file, "wb");
                    fwrite(kmerBuf[bin],kmerSize,binTally[bin],fp);
                    fclose(fp);
                    // write the tally array
                    char count_file[100];
                    sprintf(count_file,"%s/%zi-mers.%zi.%zi.count",outdir,k,bin,batches[bin]);
                    fp = fopen(count_file, "wb");
                    fwrite(tally.data(),1,distinct,fp);
                    fclose(fp);
                }
            }
            binTally[bin] = 0;
        }
    }
}

void Kmerizer::printKmer(kword_t * kmer) {
    size_t bpw = 8 * sizeof(kword_t);
    for (size_t j = 0; j < nwords; j++) {
        fprintf(stderr," ");
        for (size_t b=0;b<bpw;b++) {
            fprintf(stderr,"%d",((*(kmer+j) & (1ULL << (bpw-b-1)))) ? 1 : 0);
        }
    }
    fprintf(stderr,"\n");
}

void
Kmerizer::bitSlice(
    kword_t *    kmers,
    const size_t n,
    BitVector **      kmer_slices,
    size_t       nbits)
{
    // initialize WAH compressed bitvectors
    for (size_t i = 0; i < nbits; i++) {
        kmer_slices[i] = new BitVector(true);
    }
    bitset<256> bbit; // for keeping track of the set bits
    size_t boff[nbits]; // beginning of the current run of 1's or 0's
    memset(boff,0,nbits*sizeof(size_t)); // initialize to 0
    unsigned int bpw = 8*sizeof(kword_t); // bits per word

    // mark the set bits in the first kmer
    if (nwords>1) {
        for (size_t w=0;w<nwords;w++) {
            unsigned int count = popCount(kmers[w]);
            for (unsigned int r=1; r<=count; r++)
                bbit.set(selectBit(kmers[w],r) + w*bpw,1);
        }
        for (size_t i=1;i<n;i++) {
            // when n is large the number of different bits between kmer i and kmer i-1 is small
            // so use xor, popCount, and selectBit to identify the changed bit positions
            kword_t * kmer = kmers + i * nwords;
            kword_t * prev = kmer - nwords;
            for (size_t w=0;w<nwords;w++) {
                kword_t x = kmer[w] ^ prev[w];
                unsigned int count = popCount(x);
                for (unsigned int r = 1; r<=count; r++) {
                    unsigned int b = selectBit(x,r) + w*bpw;
                    kmer_slices[b]->appendFill(bbit.test(b),i-boff[b]);
                    bbit.flip(b);
                    boff[b] = i;
                }
            }
        }
    }
    else { // stripped down version for nwords=1
        unsigned int count = popCount(kmers[0]);
        for (unsigned int r=1; r<=count; r++)
            bbit.set(selectBit(kmers[0],r),1);
        for (size_t i=1;i<n;i++) {
            // when n is large the number of different bits between kmer i and kmer i-1 is small
            // so use xor, popCount, and selectBit to identify the changed bit positions
            kword_t x = kmers[i] ^ kmers[i-1];
            unsigned int count = popCount(x);
            for (unsigned int r = 1; r<=count; r++) {
                unsigned int b = selectBit(x,r);
                kmer_slices[b]->appendFill(bbit.test(b),i-boff[b]);
                bbit.flip(b);
                boff[b] = i;
            }
        }
    }
    // append the last runs
    for (size_t b=0; b<nbits; b++)
        if(boff[b] < n-1)
            kmer_slices[b]->appendFill(bbit.test(b),n-boff[b]);
}

void Kmerizer::doWriteBatch(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to; bin++) {
        // convert kmerBuf[bin] into a bit-sliced bitmap index
 //       fprintf(stderr,"doWriteBatch() binTally[%zi]: %u\n",bin,binTally[bin]);
        const size_t nbits = 8*kmerSize;
        BitVector* kmer_slices[nbits];
        // bitSlice(kmerBuf[bin],binTally[bin],kmer_slices,nbits);
        
        FILE *fp;
        if (0) {
            // write the kmer_slices to a file
            char kmer_file[100];
            sprintf(kmer_file,"%s/%zi-mers.%zi.%zi",outdir,k,bin,batches[bin]);
            fp = fopen(kmer_file, "wb");
            size_t n_slices = 8*sizeof(kword_t);
            fwrite(&n_slices,sizeof(size_t),1,fp);
            for (size_t b=0;b<n_slices;b++) {
                uint32_t c = kmer_slices[b]->cnt();
                fwrite(&c,sizeof(uint32_t),1,fp);
            }
            // fprintf(stderr,"binTally[%zi]: %u\n",bin,binTally[bin]);
            for (size_t b=0;b<n_slices;b++) {
                uint32_t *buf;
                size_t bytes = kmer_slices[b]->dump(&buf);
                // fprintf(stderr,"%zi: %zi,",b,bytes);
                fwrite(&bytes,sizeof(size_t),1,fp);
                fwrite(buf,1,bytes,fp);
                free(buf);
                delete kmer_slices[b];
            }
//          fprintf(stderr,"\n");
//            fprintf(stderr,"Filepointer is now at %d.\n", ftell(fp));
            fclose(fp);
        }
        // open output file for counts
        char counts_file[100];
        sprintf(counts_file,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,batches[bin]);
        fp = fopen(counts_file, "wb");

        // first write the number of distinct values
        // then write the distinct values
        // for each distinct value, write out the BitVector (size,count,rle,words.size(),words)
        size_t n_distinct = kmerFreq[bin].size();
        fwrite(&n_distinct,sizeof(size_t),1,fp);
        fwrite(kmerFreq[bin].data(),sizeof(uint32_t),n_distinct,fp);
        for (size_t i=0;i<n_distinct;i++) {
            uint32_t *buf;
            size_t bytes = counts[bin][i]->dump(&buf);
            // how many bytes are we writing
            fwrite(&bytes,sizeof(size_t),1,fp);
            // write them
            fwrite(buf,1,bytes,fp);
            free(buf);
            delete counts[bin][i]; // BitVector destructor

        }
        fclose(fp);
    }
}

// when reading counts, values is the array of kmerFrequencies
// when reading kmers, values is the number of set bits in each bit slice
void Kmerizer::readBitmap(const char* idxfile, vector<uint32_t> &values, vector<BitVector*> &index) {
    // open idxfile
    FILE *fp;
    fp = fopen(idxfile,"rb");
    // fread to repopulate values and bvecs
    size_t n_distinct;
    fread(&n_distinct,sizeof(size_t),1,fp);
    values.resize(n_distinct);
    fread(values.data(),sizeof(uint32_t),n_distinct,fp); // assumes 32-bit words in bvec
    index.resize(n_distinct);
    for (size_t i=0;i<n_distinct;i++) {
        size_t bytes;
        fread(&bytes,sizeof(size_t),1,fp);
        size_t words = bytes/sizeof(uint32_t);
        uint32_t buf[words];
        fread(buf,1,bytes,fp);
        index[i] = new BitVector(buf);
    }
    fclose(fp);
}

void Kmerizer::doLoadIndex(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to;bin++) {
        char fname[100];
        sprintf(fname,"%s/%zi-mers.%zi.idx",outdir,k,bin);
        readBitmap(fname,kmerFreq[bin],counts[bin]);
    }
}

void Kmerizer::doMergeBatches(const size_t from, const size_t to) {
    for (size_t bin=from; bin<to; bin++) {
        //fprintf(stderr,"doMergeBatches() bin: %zi, batches: %zi\n",bin,batches[bin]);
        if (lutTally[bin] > 1 || kmerLutV[0]>0) {
            
            char kmer_file[100];
            sprintf(kmer_file,"%s/%zi-mers.%zi.0.raw",outdir,k,bin);
            FILE *fp;
            fp = fopen(kmer_file, "wb");
            fwrite(kmerLutK[bin],kmerSize,lutTally[bin],fp);
            fclose(fp);
            // write the tally array
            char count_file[100];
            sprintf(count_file,"%s/%zi-mers.%zi.0.count",outdir,k,bin);
            fp = fopen(count_file, "wb");
            fwrite(kmerLutV[bin],sizeof(uint32_t),lutTally[bin],fp);
            fclose(fp);
        }

/*
        char fname[100];
        // initialize the BitSlicedIndex for the common kmers in kmerLutK
        sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);

        BitSlicedIndex *mergedKmers = new BitSlicedIndex(nwords,fname);
        
        if (batches[bin] == 0) {
            // only have to deal with the common kmers from kmerLutK
            for (int i=0;i<lutTally[bin];i++) {
                mergedKmers->append(kmerLutK + i*nwords);
            }
            // create a range encoded bitmap index for the common kmer counts in kmerLutV
            sprintf(fname,"%s/%zi-mers.%zi.idx",outdir,k,bin);
            RangeEncodedIndex *kmerCounts = new RangeEncodedIndex(kmerLutV,fname);
            kmerCounts->saveIndex();
        }
        else {
            // need to merge in-memory lookup table with serialized batches
            // open each file of raw kmers
            // try mmap-ing so we can treat them like arrays
            // the counts were serialized as an array of uint8_t
            uint8_t *kmerCounts [batches[bin]];
            kword_t *kmerValues [batches[bin]]; // mmap'ed so memory is managed by kernel
            for (int i=0;i<batches[bin];i++) {
                // load the counts from batch i
                sprintf(fname,"%s/%zi-mers.%zi.%i.count",outdir,k,bin,i);
                kmerCounts[i];

                // mmap the raw kmers file from batch i
                sprintf(fname,"%s/%zi-mers.%zi.%i.raw",outdir,k,bin,i);
                kmerValues[i];
            }
            vector<uint32_t> mergedCounts;
            // setup a heap with the min value from each batch
            // initialize super fast direct hash table
            // initialize the cutoff
            while (heap not empty) {
                // extract the min from the heap
                // if < cutoff, take any other values from the corresponding batch below cutoff
                // add each to the counting hash
                // insert the next kmer into the heap
                // once you hit a key that is too big, stop
                // iterate over the set bits and output keys to mergedKmers and counts to mergedCounts
                // zero the hash while processing
                mergedKmers->append();
                mergedCounts.push_back();
            }

            // create the range encoded index for the mergedCounts
            sprintf(fname,"%s/%zi-mers.%zi.idx",outdir,k,bin);
            RangeEncodedIndex *kmerCounts = new RangeEncodedIndex(mergedCounts,fname);
            kmerCounts->saveIndex();
        }
        mergedKmers->saveIndex();
*/
    }
}
/*
        // read the counts and kmers for each batch
        vector<BitVector*>    batch_counts[batches[bin]];
        vector<uint32_t> batch_values[batches[bin]];
        kword_t  kmers[batches[bin]*nwords]; // next kmer in each active batch
        uint32_t btally[batches[bin]]; // frequency of next kmer in each batch
        size_t   offset[batches[bin]]; // keep track of position in each batch

        for (size_t i=0;i<batches[bin];i++) {
            char fname[100];
            sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i+1);
            readBitmap(fname,batch_values[i],batch_counts[i]);
        }

        memset(offset, 0, sizeof(size_t)*batches[bin]);
        uint32_t todo = batches[bin]; // number of batches to process
        vector<uint32_t> tally;

        const size_t nbits = 8 * kmerSize;
        BitVector* merged_slices[nbits];
        bitset<256> bbit;
        unsigned int boff[nbits];
        unsigned int n=0;
        memset(boff, 0, nbits * sizeof(size_t));
        const unsigned int bpw = 8 * sizeof(kword_t); // bits per word

        for (size_t b=0;b<nbits;b++)
            merged_slices[b] = new BitVector(true);

        // read the first kmer and counts from each batch
        for (size_t i=0;i<batches[bin];i++) {
            size_t err =
                pos2kmer(offset[i], kmers + i * nwords, batch_slices[i]);
            if (err != 0) {
                todo--;
                btally[i]=0;
            }
            else {
                btally[i] = pos2value(offset[i]++,batch_values[i],batch_counts[i]);
            }
        }
        // choose min
        size_t mindex = findMin(kmers,btally);
        kword_t distinct[nwords];
        memcpy(distinct,kmers + mindex*nwords, kmerSize);
        tally.push_back(btally[mindex]);
        // mark the set bits in the first distinct kmer
        for (size_t w=0;w<nwords;w++) {
            unsigned int count = popCount(distinct[w]);
            for (unsigned int r=1; r<=count; r++)
                bbit.set(selectBit(distinct[w],r)+w*bpw,1);
        }
        // replace min
        // iterate until there's nothing left to do
    
        while(todo>0) {
            size_t err = pos2kmer(offset[mindex],kmers + mindex*nwords, batch_slices[mindex]);
            if (err != 0) {
                todo--;
                btally[mindex]=0;
                if (todo==0) break;
            }
            else
                btally[mindex] = pos2value(offset[mindex]++,batch_values[mindex],batch_counts[mindex]);
            mindex = findMin(kmers,btally);
            // compare to distinct
            if (kmercmp(kmers + mindex*nwords, distinct, nwords) == 0) // same kmer
                tally.back() += btally[mindex];
            else { // find the changed bits
                tally.push_back(btally[mindex]);
                n++;
                kword_t *minkmer = kmers + mindex*nwords;
                for (size_t w=0;w<nwords;w++) {
                    kword_t x = distinct[w] ^ minkmer[w];
                    unsigned int count = popCount(x);
                    for (unsigned int r = 1; r<=count; r++) {
                        unsigned int b = selectBit(x,r) + w*bpw;
                        merged_slices[b]->appendFill(bbit.test(b),n-boff[b]);
                        bbit.flip(b);
                        boff[b]=n;
                    }
                    distinct[w] = minkmer[w];
                }
            }
        }
        // finish the bitvectors
        n++;
        for (size_t b=0; b<nbits; b++)
            if(boff[b]<n-1)
                merged_slices[b]->appendFill(bbit.test(b),n-boff[b]);
    
        rangeIndex(tally,kmerFreq[bin],counts[bin]);

        // clean up the bitvector pointers
        for(size_t i=0;i<batches[bin];i++)
            for (size_t b=0;b<nbits;b++) {
                delete batch_slices[i][b];
                delete batch_counts[i][b];
            }

        // open an output file for the merged distinct kmers
        FILE *ofp;
        char fname[100];
        sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
        ofp = fopen(fname,"wb");

        fwrite(&nbits,sizeof(size_t),1,ofp);
        for (size_t b=0;b<nbits;b++) {
            uint32_t c = merged_slices[b]->cnt();
            fwrite(&c,sizeof(uint32_t),1,ofp);
        }
        for (size_t b=0;b<nbits;b++) {
            uint32_t *buf;
            size_t bytes = merged_slices[b]->dump(&buf);
            // fprintf(stderr,"%zi: %zi,",b,bytes);
            fwrite(&bytes,sizeof(size_t),1,ofp);
            fwrite(buf,1,bytes,ofp);
            free(buf);
            delete merged_slices[b];
        }

        fclose(ofp);
    
        char counts_file[100];
        sprintf(counts_file,"%s/%zi-mers.%zi.idx",outdir,k,bin);
        ofp = fopen(counts_file, "wb");
        // first write the number of distinct values
        // then write the distinct values
        // for each distinct value, write out the BitVector (size,count,rle,words.size(),words)
        int n_distinct = kmerFreq[bin].size();
        fwrite(&n_distinct,sizeof(int),1,ofp);
        fwrite(kmerFreq[bin].data(),sizeof(uint32_t),n_distinct,ofp);
        for (size_t i=0;i<n_distinct;i++) {
            uint32_t *buf;
            size_t bytes = counts[bin][i]->dump(&buf);
            fwrite(&bytes,sizeof(size_t),1,ofp);
            fwrite(counts[bin][i]->getWords().data(),1,counts[bin][i]->bytes(),ofp);
            delete counts[bin][i];
        }
        fclose(ofp);
    
        // delete input files
        for (size_t i=0;i<batches[bin];i++) {
            sprintf(fname,"%s/%zi-mers.%zi.%zi",outdir,k,bin,i+1);
            if (remove(fname) != 0) perror("error deleting file");
            sprintf(fname,"%s/%zi-mers.%zi.%zi.idx",outdir,k,bin,i+1);
            if (remove(fname) != 0) perror("error deleting file");
        }
    }
}
*/
size_t Kmerizer::pos2kmer(size_t pos, kword_t *kmer, vector<BitVector*> &index) {
    if (pos >= index[0]->getSize()) return 1;
    size_t bpw = 8*sizeof(kword_t);
    // need to zero the kmer first
    memset(kmer,0,kmerSize);
    for(size_t b=0;b<8*kmerSize;b++)
        if (index[b]->find(pos))
            kmer[b/bpw] |= 1ULL << (bpw - b%bpw - 1);
    return 0;
}

uint32_t Kmerizer::pos2value(
    size_t pos,
    vector<uint32_t> &values,
    vector<BitVector*> &index)
{
    // lookup the value in the pos bit
    // find the first BitVector where this bit is set
    for(size_t i=0;i<values.size();i++)
        if (index[i]->find(pos))
            return values[i];
    return 0;
}


inline size_t
Kmerizer::findMin(const kword_t * kmers, const uint32_t * kcounts) {
    size_t mindex = 0;
    // while (kcounts[mindex] == 0) mindex++;
    // for (size_t i=mindex+1;i<batches[bin];i++)
    //     if (kcounts[i] != 0)
    //         if (kmercmp(kmers + i*nwords,kmers + mindex*nwords,nwords) < 0)
    //             mindex = i;

    return mindex;
}

void Kmerizer::vecHist(vector<uint32_t> &vec, vector<uint32_t> &values, vector<uint32_t> &frequency) {
//    fprintf(stderr,"vecHist()");
    timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);

    vector<uint32_t> tmp = vec;
    sort(tmp.begin(),tmp.end());
    values.clear();
    frequency.clear();
    int offset=0;
    values.push_back(tmp[0]);
    frequency.push_back(1);
    for(int i=1;i<tmp.size();i++) {
        if (tmp[i] == values[offset]) {
            frequency[offset]++;
        }
        else {
            offset++;
            frequency.push_back(1);
            values.push_back(tmp[i]);
        }
    }
    gettimeofday(&t2, NULL);
    elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
//    fprintf(stderr," took %f seconds\n",elapsedTime);
}

void print_binary(uint64_t w) {
    for (int b = 63; b >= 0; b--)
        fprintf(stderr,"%d", w & (1ULL << b) ? 1 : 0);
    fprintf(stderr,"\n");
}

void Kmerizer::bitHist(vector<uint32_t> &vec, vector<uint32_t> &values, vector<uint32_t> &frequency) {
    // fprintf(stderr,"bitHist()");
    // timeval t1, t2;
    // double elapsedTime;
    // gettimeofday(&t1, NULL);
    values.clear();
    kword_t mybits [256]; // actually 64*256 = 16384 bits 2^14
    memset(mybits,0,256*8);
    vector<uint32_t> overflow;
    vector<uint32_t>::iterator it;
    for (it = vec.begin(); it < vec.end(); ++it)
        if (*it < 16384UL) // first 8 bits choose the word, last 6 bits choose the bit
            mybits[*it >> 6] |= (1ULL << (63 - (*it & 63UL)));
        else
            overflow.push_back(*it);
    // find the set bits
    // save the popCount for use later
    unsigned int totalpc[256];
    for (unsigned int i=0; i<256; i++) {
        unsigned int pc = popCount(mybits[i]);
        for (unsigned int r=1; r<=pc; r++) {
            unsigned int b = selectBit(mybits[i], r);
            values.push_back((i << 6) | b);
        }
        if (i>0)
            totalpc[i] = totalpc[i-1] + pc;
        else
            totalpc[i] = pc;
    }
    unsigned int smallCount = values.size();
    if (smallCount != totalpc[255]) {
        fprintf(stderr,"total popCount (%u) doesn't match size of values vector (%zi)\n",totalpc[255],values.size());
    }
//    fprintf(stderr,"smallCount: %i, overflow.size(): %zi\n",smallCount,overflow.size());
    // uniqify the overflow
    frequency.clear();
    frequency.resize(values.size());
    if (overflow.size() > 0) {
        sort(overflow.begin(),overflow.end());
        frequency.push_back(1);
        values.push_back(overflow[0]);
        for(int i=1;i<overflow.size();i++) {
            if (overflow[i] == values.back())
                frequency.back()++;
            else {
                values.push_back(overflow[i]);
                frequency.push_back(1);
            }
        }
    }
    // gettimeofday(&t2, NULL);
    // elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    // fprintf(stderr," uniqifying took %f seconds. ",elapsedTime);

    // count the frequency
    for (it = vec.begin(); it < vec.end();++it) {
        if (*it < 16384) {
            unsigned int bin = *it >> 6;
            int offset = 0;
            if(bin>0)
                offset = totalpc[bin-1];
            // count set bits in mybits[bin] before this bit
            if (*it & 63UL) {
                // print_binary(mybits[bin]);
                // print_binary(0xFFFFFFFFFFFFFFFFULL << 64 - (*it & 63UL));
                // print_binary(mybits[bin] & (0xFFFFFFFFFFFFFFFFULL << 64 - (*it & 63UL)));
                offset += popCount(mybits[bin] & (0xFFFFFFFFFFFFFFFFULL << 64 - (*it & 63UL)));
            }
            // fprintf(stderr,"%u bin: %i, : %i\n",*it,bin,offset);
            frequency[offset]++;
        }
    }
    // gettimeofday(&t2, NULL);
    // elapsedTime = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    // fprintf(stderr," counting took %f seconds\n",elapsedTime);
}
    
// for each distinct value in the vec create a bitvector
// indexing the positions in the vec holding a value <= v

// vector<uint32_t> values would be more compact as a bitvector.

// test a diy bitset class that uses the fast popcount,rank, and select functions
// there could be dense/sparse regions in the bitset described by an aux array
// for this purpose you would probably have one small bitset followed by a vector of bigger values
void Kmerizer::rangeIndex(vector<uint32_t> &vec, vector<uint32_t> &values, vector<BitVector*> &index) {
    // find the distinct values in vec
    values.clear();
    // use a bitset to mark the distinct values (from 0-255)
    // and another vector for overflow (which we'll sort/uniq later)
    bitset<256> mybits;
    vector<uint32_t> overflow;
    vector<uint32_t>::iterator it;
    for (it = vec.begin(); it < vec.end(); ++it) {
        if (*it >= 256) 
            overflow.push_back(*it);
        else
            mybits[*it]=1;
    }
    // populate values vector with set bits in mybits
    for (uint32_t i=0;i<256;i++)
        if (mybits.test(i))
            values.push_back(i);
    
    if (overflow.size() > 0) {
        sort(overflow.begin(),overflow.end());
        it = unique(overflow.begin(),overflow.end());
        values.insert(values.end(),overflow.begin(),it);
    }
//    fprintf(stderr,"rangeIndex() %zi distinct values\n",values.size());
    // setup a vector for each range
    vector<uint32_t> vrange[values.size()];
    // iterate over the vec and push the offset onto each range
    for (size_t i=0;i<vec.size();i++) {
        it = lower_bound(values.begin(),values.end(),vec[i]);
        for (size_t j=0;j<=(it - values.begin());j++)
            vrange[j].push_back(i);
    }
    // create a bitvector for each range
    index.resize(values.size());
    for (size_t i=0;i<values.size();i++)
        index[i] = new BitVector(vrange[i]);
}

// need an iterator that works with a BitVector mask BitVector->next_one(); - returns position of set bit
// or number of bits (invalid response)

// We can do this in parallel if we figure out in advance the offset within the output file for each bin.
// A frequency histogram can get us the number of 1 digit counts, 2 digit counts, 3 digit counts, etc.
// So, if there are mask->cnt() set bits, we need n*(k+2) + n1 + 2n2 + 3n3 etc bytes. Each kmer takes k bytes, and you have a tab and a newline character on each line. n1 = 1 digit counts, n2 = 2 digit counts, etc...
void Kmerizer::dump(char *fname) {
    // create a bitvector of all 1's for each bin
    // and call dump(FILE *fp, BitVector *mask)
    BitVector *mask[NBINS];
    for (size_t i=0;i<NBINS;i++) {
        mask[i] = new BitVector(true); // first create an empty bitvector
        mask[i]->appendFill(true,counts[i][0]->getSize()); // then make it all 1's
    }
    sdump(fname,mask);
}

// find kmers with frequencies in the given range[min,max]
// when max < min, ignore max and find kmers with frequencies in the range[min,infinity]
void Kmerizer::filter(uint32_t min, uint32_t max, BitVector **mask) {
    // do range queries over each bin to build up an array of BitVector masks.
    boost::thread_group tg;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        tg.create_thread(boost::bind(&Kmerizer::doFilter, this, i, j, min, max, mask));
    }
    tg.join_all();
}


void Kmerizer::pdump(char *fname, BitVector **mask) {
    // determine the uncompressed bin sizes in advance so we can doDump in parallel
    // binsize[bin] = mask[bin]->cnt() * (k+2) + d*number of d digit counts
    uint32_t binsize[NBINS];
    uint32_t total_size=0;
    for (size_t bin=0;bin<NBINS;bin++) {
        binsize[bin] = mask[bin]->cnt() * (k+2); // k-mer\t\n
        // d digit counts
        uint32_t d=1;
        vector<uint32_t>::iterator it = lower_bound(kmerFreq[bin].begin(),kmerFreq[bin].end(),10);
        while(it != kmerFreq[bin].end()) {
            if (it > kmerFreq[bin].begin()) {
                it--;
                BitVector *dmask = *counts[bin][*it] & *mask[bin];
                binsize[bin] += d*dmask->cnt();
            }
            d++;
            it = lower_bound(it,kmerFreq[bin].end(),10^d);
        }
        total_size += binsize[bin];
    }
    // use memory mapped IO to do this
    // maybe open the output file, resize it, and then mmap it to some buffer
    // You can't do this with stdout. is that a problem?
    char *buff;
//  mmap(buff, total_size);
    boost::thread_group tg;
    uint32_t offset=0;
    for (size_t i = 0; i < NBINS; i += threadBins) {
        size_t j = (i + threadBins > NBINS) ? NBINS : i + threadBins;
        // pass a pointer to the right place in the mmap'd output file
        tg.create_thread(boost::bind(&Kmerizer::doPdump, this, i, j, buff + offset, mask));
        for (size_t bin=i;bin<j;bin++)
            offset += binsize[bin];
    }
    tg.join_all();
    // close the output file
}

void Kmerizer::sdump(char *fname, BitVector **mask) {
    // open output file
    FILE *fp;
    fp = fopen(fname, "w");
    size_t bpw = 8*sizeof(kword_t); // bits per word
    char kstr[k+1]; // unpack each kmer into this char array.
    for (size_t bin=0;bin<NBINS;bin++) {
        // check if bit slices have been loaded
        if (slices[bin].empty()) {
            char bfname[100];
            sprintf(bfname,"%s/%zi-mers.%zi",outdir,k,bin);
            vector<uint32_t> tmp;
            readBitmap(bfname,tmp,slices[bin]);
        }
        // iterate over the set bits in mask[bin]
        uint32_t bv_len = mask[bin]->getSize();
        uint32_t next1 = mask[bin]->nextOne(0);
        while (next1 < bv_len) {
            kword_t kmer[nwords];
            for (size_t w=0;w<nwords;w++) {
                kmer[w]=0;
                for (size_t b=0;b<bpw;b++) {
                    if(slices[bin][w*bpw + b]->find(next1)) {
                        kmer[w] |= 1 << bpw-b-1;
                    }
                }
            }
            // output the kmer and count
            unpack(kmer, kstr);
            fprintf(fp,"%s\t%u\n",kstr,frequency(bin,next1));
            uint32_t next1 = mask[bin]->nextOne(next1+1);
        }
    }
    fclose(fp);
}

void
Kmerizer::doPdump(const size_t from,
                 const size_t to,
                 char *       buff,
                 BitVector ** mask)
{
    size_t bpw = 8 * sizeof(kword_t); // bits per word
    char kstr[k+1]; // unpack each kmer into this char array.
    for (size_t bin = from; bin < to; bin++) {
        // check if we've loaded the bit slices for this bin
        if (slices[bin].empty()) {
            // load the kmer slices for this bin
            char fname[100];
            sprintf(fname,"%s/%zi-mers.%zi",outdir,k,bin);
            vector<uint32_t> junk;
            readBitmap(fname,junk,slices[bin]);
        }
        // iterate over the set bits in mask[bin]
        uint32_t bv_len = mask[bin]->getSize();
        uint32_t next1 = mask[bin]->nextOne(0);
        while (next1 < bv_len) {
            kword_t kmer[nwords];
            for (size_t w=0;w<nwords;w++) {
                kmer[w]=0;
                for (size_t b=0;b<bpw;b++) {
                    if(slices[bin][w*bpw + b]->find(next1)) {
                        kmer[w] |= 1 << bpw-b-1;
                    }
                }
            }
            // output the kmer and count
            unpack(kmer, kstr);
            // instead of fprintf, use memcpy() to write kstr and frequency to the mmap'd output file
//          fprintf(fp,"%s\t%u\n",kstr,frequency(bin,next1));
            uint32_t next1 = mask[bin]->nextOne(next1+1);
        }
    }
}

// range encoding: bitvector x marks kmers that occur <= x times
// therefore counts[bin].back() marks all rows
void
Kmerizer::doFilter(const size_t from,
                   const size_t to,
                   uint32_t min,
                   uint32_t max,
                   BitVector **mask)
{
    for (size_t bin = from; bin < to; bin++) {
        if (max > 0) {
            vector<uint32_t>::iterator it =
                lower_bound(kmerFreq[bin].begin(),kmerFreq[bin].end(),max);
            if (it == kmerFreq[bin].begin()) {
                // create a null mask and set the length
                mask[bin] = new BitVector(true);
                mask[bin]->appendFill(false,counts[bin][0]->getSize());
            }
            else {
                if (*it < max) it--;
                mask[bin]->copy(counts[bin][it - kmerFreq[bin].begin()]);
            }
        }
        else
            mask[bin]->copy(counts[bin].back()); // this is how we deal with max <= 0 (no max)
        if (mask[bin]->cnt() > 0 && min>1) { // nand
            // need to find the first BitVector thats < min
            vector<uint32_t>::iterator it = lower_bound(kmerFreq[bin].begin(),kmerFreq[bin].end(),min);
            if (it > kmerFreq[bin].begin()) {
                it--;
                *mask[bin] &= *(counts[bin][it - kmerFreq[bin].begin()]->copyflip());
            }
        }
    }
}
