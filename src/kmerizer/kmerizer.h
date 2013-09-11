#ifndef SNAPDRAGON_KMERIZER_H
#define SNAPDRAGON_KMERIZER_H

#define DEBUG true
#define NBINS 256
#define READ 1
#define QUERY 2
#define CANONICAL 'C'

#include <vector>
#include "../bvec/bvec.h"

typedef uint64_t kword_t;
using namespace std;

class Kmerizer {
    size_t  k;
    kword_t kmask;
    size_t  shiftlastby;
    size_t  rshift, lshift;
    size_t  nwords;
    size_t  kmerSize; // in bytes
    size_t  threads;
    size_t  threadBins;
    size_t  batches[NBINS];
    char    mode;
    char    state;
    char *  outdir;

    // unsorted padded packed kmers
    kword_t *             kmerBuf[NBINS];
    
    // number of kmers in each bin
    uint32_t              binTally[NBINS];
    uint32_t              rareKmers[NBINS]; // number of rare kmers in each bin (not in kmerLut)

    // capacity of each bin
    uint32_t              binCapacity[NBINS];
    uint32_t              totalCapacity[NBINS];
    
    // common kmers lookup table: implemented as a pair of arrays
    kword_t *             kmerLutK[NBINS]; // kmer
    uint32_t *            kmerLutV[NBINS]; // frequency
    uint32_t              lutTally[NBINS]; // number of entries

    // data structures for querying
    // sorted distinct kmer frequencies
    vector<uint32_t>      kmerFreq[NBINS];

    // bitmap index of frequency counts
    vector<BitVector*>    counts[NBINS];

    // bit sliced bitmap self index of sorted kmers
    vector<BitVector*>    slices[NBINS];

public:
    // constructor
    Kmerizer(const size_t k,
             const size_t threads,
             const char * outdir,
             const char   mode);
    
    ~Kmerizer() {};

    // allocate memory for internal data structures
    int allocate(const size_t maximem);

    // extract (canonicalized) kmers from the sequence
    void addSequence(const char* seq,const int length);

    // write distinct kmers and RLE counts to disk (merging multiple batches)
    void save();

    // read kmer indexes into memory
    void load();

    // output the kmer count frequency distribution
    void histogram();

    // output the frequency of the kmers in the given sequence
    uint32_t find(const char* query);

    // write the kmers to an output file as a tab delimited dump (kmer, frequency)
    void dump(char *fname);
    
    // range query - populates mask bitvectors to mark kmers in range
    void filter(uint32_t min, uint32_t max, BitVector **mask);

    void vecHist(vector<uint32_t> &vec, vector<uint32_t> &values, vector<uint32_t> &frequency);
    void bitHist(vector<uint32_t> &vec, vector<uint32_t> &values, vector<uint32_t> &frequency);

private:

    // pack nucleotides into 2 bits
    inline kword_t twoBit(const kword_t val) const;

    // shift kmer to make room for nucl
    size_t nextKmer(kword_t* kmer, size_t bin, const char nucl);

    inline void unpack(kword_t* kmer, char *seq);

    kword_t* canonicalize(kword_t *packed, kword_t *rcpack, size_t *bin) const;

    // reverse complement
    inline kword_t revcomp(const kword_t val) const;

    inline int searchLut1(kword_t *kmer, size_t bin);
    int searchLut(kword_t *kmer, size_t bin);

    // count set bits in a kmer (usually XOR of 2 kmers)
    inline unsigned int popCount(kword_t v) const;

    // returns the position of the rth set bit in v
    inline unsigned int selectBit(kword_t v, unsigned int r);

    size_t pos2kmer(size_t pos, kword_t *kmer, vector<BitVector*> &index);

    uint32_t pos2value(size_t pos, vector<uint32_t> &values,
                       vector<BitVector*> &index);

    uint32_t find(kword_t *kmer, size_t bin);

    // kmerBuf is full. uniqify and write batch to disk
    void serialize();
    
    // sort each kmerBuf, update binTally, and fill counts
    void uniqify();
    void doUnique(const size_t from, const size_t to);
    // write the counted kmers to an output file
    void writeBatch();
    void doWriteBatch(const size_t from, const size_t to);
    // merge output files
    void mergeBatches();
    void doMergeBatches(const size_t from, const size_t to);
    // read the frequency bitmap index into memory
    void doLoadIndex(const size_t from, const size_t to);
    // write kmers to an output file
    void doDump(const size_t from, const size_t to, FILE *fp, BitVector **mask);
    // range query
    void doFilter(const size_t from,
                  const size_t to,
                  uint32_t min,
                  uint32_t max,
                  BitVector **mask);
    
    void pdump(char *fname, BitVector **mask);
    void sdump(char *fname, BitVector **mask);
    void doPdump(const size_t from,
                 const size_t to,
                 char *buff,
                 BitVector **mask);
    // is this too generic to go here?
    void rangeIndex(vector<uint32_t> &vec,
                    vector<uint32_t> &values,
                    vector<BitVector*> &index);
    void readBitmap(const char* idxfile,
                    vector<uint32_t> &values,
                    vector<BitVector*> &index);
    size_t findMin(const kword_t* kmers, const uint32_t* kcounts);
    // uint32_t pos2value(size_t pos, vector<uint32_t> &values, vector<BitVector*> &index);
    // size_t pos2kmer(size_t pos, kword_t *kmer, vector<BitVector*> &index);
    // lookup the frequency of a specific kmer by position?

    uint32_t frequency(size_t bin, uint32_t pos);

    void printKmer(kword_t *kmer);
    void bitSlice(kword_t *kmers,
                  const size_t n,
                  BitVector **kmer_slices,
                  size_t nbits);

    
};

inline kword_t Kmerizer::twoBit(const kword_t val) const {
    static const kword_t table[256] =
    {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,
        0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    return table[val];
}

// quickly count the number of set bits in a k-mer
inline unsigned int Kmerizer::popCount(kword_t v) const {
    // number of 1 bits in a value between 0 and 255
    static const unsigned char t[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
    unsigned char * p = (unsigned char *) &v;
    return t[p[0]]+t[p[1]]+t[p[2]]+t[p[3]]+t[p[4]]+t[p[5]]+t[p[6]]+t[p[7]];
}

// select the bit position given the rank
inline unsigned int Kmerizer::selectBit(kword_t v, unsigned int r) {
    unsigned int s;
    kword_t a, b, c, d;
    unsigned int t;
    // a = (v & 0x5555...) + ((v >> 1) & 0x5555...);
    a =  v - ((v >> 1) & ~0UL/3);
    // b = (a & 0x3333...) + ((a >> 2) & 0x3333...);
    b = (a & ~0UL/5) + ((a >> 2) & ~0UL/5);
    // c = (b & 0x0f0f...) + ((b >> 4) & 0x0f0f...);
    c = (b + (b >> 4)) & ~0UL/0x11;
    // d = (c & 0x00ff...) + ((c >> 8) & 0x00ff...);
    d = (c + (c >> 8)) & ~0UL/0x101;
    t = (d >> 32) + (d >> 48);
    // Now do branchless select!                                                
    s  = 64;
    // if (r > t) {s -= 32; r -= t;}
    s -= ((t - r) & 256) >> 3; r -= (t & ((t - r) >> 8));
    t  = (d >> (s - 16)) & 0xff;
    // if (r > t) {s -= 16; r -= t;}
    s -= ((t - r) & 256) >> 4; r -= (t & ((t - r) >> 8));
    t  = (c >> (s - 8)) & 0xf;
    // if (r > t) {s -= 8; r -= t;}
    s -= ((t - r) & 256) >> 5; r -= (t & ((t - r) >> 8));
    t  = (b >> (s - 4)) & 0x7;
    // if (r > t) {s -= 4; r -= t;}
    s -= ((t - r) & 256) >> 6; r -= (t & ((t - r) >> 8));
    t  = (a >> (s - 2)) & 0x3;
    // if (r > t) {s -= 2; r -= t;}
    s -= ((t - r) & 256) >> 7; r -= (t & ((t - r) >> 8));
    t  = (v >> (s - 1)) & 0x1;
    // if (r > t) s--;
    s -= ((t - r) & 256) >> 8;
    s = 65 - s;
    return s-1;
}

inline void Kmerizer::unpack(kword_t* kmer, char* seq) {
    static const char table[4] = {65, 67, 71, 84};
    for(size_t i=0;i<k;i++) {
        // which word, which nucl
        size_t w = i>>5;
        if (w == nwords-1) // not a full length word
            seq[i] = table[(kmer[w] >> 2*(k-i)-2) & 3];
        else
            seq[i] = table[(kmer[w] >> (62 - 2*(i%32))) & 3];
    }
    seq[k] = '\0';
}

inline kword_t Kmerizer::revcomp(const kword_t val) const {
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

    return
        (rctable[val&0xFFUL]<<56) |
        (rctable[(val>>8)&0xFFUL]<<48) |
        (rctable[(val>>16)&0xFFUL]<<40) |
        (rctable[(val>>24)&0xFFUL]<<32) |
        (rctable[(val>>32)&0xFFUL]<<24) |
        (rctable[(val>>40)&0xFFUL]<<16) |
        (rctable[(val>>48)&0xFFUL]<<8) |
        (rctable[(val>>56)&0xFFUL]);
}


// diy binary search
inline int Kmerizer::searchLut1(kword_t *kmer, size_t bin) {
    kword_t* A = kmerLutK[bin];
    int imin = 0;
    int imax = lutTally[bin];
    while (imin < imax) {
        int imid = (imin+imax)>>1;
//        int cmp = kmercmp(A + nwords*imid, kmer, nwords);
        kword_t *mid = A + nwords*imid;
        if (*mid < *kmer)
            imin = imid + 1;
        else if (*mid > *kmer)
            imax = imid;
        else
            return imid;
    }
    return -1;
}


#endif // #ifndef SNAPDRAGON_KMERIZER_H
