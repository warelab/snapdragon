#ifndef KMERIZER_H
#define KMERIZER_H

#define DEBUG false
#define NBINS 256
#define CANONICAL 'C'
#define BOTH 'B'
#define READING 1
#define QUERY 2

#include <vector>
#include "../bvec/bvec.h"

typedef uint64_t kword_t;
using namespace std;

class kmerizer {
    size_t k;
    kword_t kmask;
    size_t shiftlastby;
    size_t rshift, lshift;
    size_t nwords;
    size_t kmer_size; // in bytes
    size_t threads;
    size_t thread_bins;
    size_t batches;
    size_t max_kmers_per_bin;
    char mode;
    char state;
    char *outdir;

    kword_t* kmer_buf[NBINS]; // raw unsorted padded kmers, or sort|uniq'ed kmers
    uint32_t bin_tally[NBINS]; // number of kmers in each bin (or number of distinct kmers)
    vector<uint32_t> kmer_freq [NBINS]; // sorted distinct kmer frequencies
    vector<bvec*> counts [NBINS]; // bitmap index of frequency counts
    vector<bvec*> slices [NBINS]; // bitmap self index of kmers

public:
    kmerizer(const size_t _k, const size_t _threads, char* _outdir, const char _mode);
    int allocate(const size_t maximem); // allocates memory for each kmer_buf
    void addSequence(const char* seq,const int length); // extract (canonicalized) kmers from the sequence
    void save(); // writes distinct kmers and rle counts to disk (merging multiple batches)
    void load(); // reads kmer indexes into memory
    void histogram(); // output the kmer count frequency distribution
    uint32_t find(char* query);
    void dump(char *fname);
    void dump(char *fname,bvec **mask);
    void filter(uint32_t min, uint32_t max, bvec **mask);
    ~kmerizer() {};

private:
    inline kword_t twobit(const kword_t val) const; // pack nucleotides into 2 bits
    void next_kmer(kword_t* kmer, const char nucl); // shift kmer to make room for nucl
    inline void unpack(kword_t* kmer, char *seq);
    inline kword_t revcomp(const kword_t val) const; // reverse complement
    inline uint8_t hashkmer(const kword_t *kmer, const uint8_t seed) const; // to select a bin
    inline unsigned int popcount(kword_t v) const; // count set bits in a kmer (actually XOR of 2 kmers)
    inline unsigned int selectbit(kword_t v, unsigned int r); // returns the position of the rth set bit in v
    inline size_t pos2kmer(size_t pos, kword_t *kmer, vector<bvec*> &index);
    inline uint32_t pos2value(size_t pos, vector<uint32_t> &values, vector<bvec*> &index);
    
    kword_t* canonicalize(kword_t *packed, kword_t *rcpack) const;
    uint32_t find(kword_t *kmer, size_t bin);
    void serialize(); // kmer_buf is full. uniqify and write batch to disk
    void uniqify(); // qsort each kmer_buf, update bin_tally, and fill counts
    void do_unique(const size_t from, const size_t to); // for parallelization
    void writeBatch();
    void do_writeBatch(const size_t from, const size_t to); // for parallelization
    void mergeBatches();
    void do_mergeBatches(const size_t from, const size_t to);
    void do_loadIndex(const size_t from, const size_t to);
    void do_dump(const size_t from, const size_t to, FILE *fp, bvec **mask);
    void do_filter(const size_t from, const size_t to, uint32_t min, uint32_t max, bvec **mask);

    // is this too generic to go here?
    void range_index(vector<uint32_t> &vec, vector<uint32_t> &values, vector<bvec*> &index);
    void read_bitmap(const char* idxfile, vector<uint32_t> &values, vector<bvec*> &index);
    size_t find_min(const kword_t* kmers, const uint32_t* kcounts);
    // uint32_t pos2value(size_t pos, vector<uint32_t> &values, vector<bvec*> &index);
    // size_t pos2kmer(size_t pos, kword_t *kmer, vector<bvec*> &index);
    void print_kmer(kword_t *kmer);
    void bit_slice(kword_t *kmers, const size_t n, bvec **kmer_slices, size_t nbits);

};

inline kword_t kmerizer::twobit(const kword_t val) const {
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
inline unsigned int kmerizer::popcount(kword_t v) const {
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
inline unsigned int kmerizer::selectbit(kword_t v, unsigned int r) {
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

inline void kmerizer::unpack(kword_t* kmer, char* seq) {
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

inline kword_t kmerizer::revcomp(const kword_t val) const {
    // bitwise reverse complement of values from 0 to 255
    static const kword_t table[256] =
    {
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
        (table[val&0xFFUL]<<56) |
        (table[(val>>8)&0xFFUL]<<48) |
        (table[(val>>16)&0xFFUL]<<40) |
        (table[(val>>24)&0xFFUL]<<32) |
        (table[(val>>32)&0xFFUL]<<24) |
        (table[(val>>40)&0xFFUL]<<16) |
        (table[(val>>48)&0xFFUL]<<8) |
        (table[(val>>56)&0xFFUL]);
}

inline uint8_t kmerizer::hashkmer(const kword_t *kmer, const uint8_t seed) const {
    static const uint8_t Rand8[256] =
    {
        105,193,195, 26,208, 80, 38,156,128,  2,101,205, 75,116,139, 61,
        197,120,244, 51,185,132, 55,150,177,241,103,196, 13,237,136, 24,
        211, 56,207,  9, 30,145, 18,167,108, 32,106,151,  3, 54,248, 65,
         17,198, 85, 95, 29, 83,  4,206,188,186,107,255,129, 35,142, 91,
        203,158, 74,138,162,135,102,114, 81,170, 10, 19,215, 57,214, 70,
         37,163,231,227,152, 14, 40, 84, 68,252,  0, 66,121,127,223, 78,
        201, 34,225,240,124,191, 12, 42, 92,209,  1, 31,130,100, 28,224,
        161,249,110, 77, 87,144,181, 21, 86, 58,174,113,194,147,242, 50,
         59, 48,250, 88,184,245, 45, 44,148, 73,154,230,149, 89,118,119,
         79,229,239,117,189,179,254,155, 20,176,157,212, 36,123,234, 46,
        159,202,171, 67, 93, 62, 47,164,247, 15,137,235,216,160,200,133,
        140,172,192,221,131,111,218,210,153,219, 41, 72, 63,183, 39,  8,
         98,168, 52,213,175,134,115, 90, 82, 16,226,220,251, 69,243,233,
        180,  5, 99, 60,  7,204,112,253,182,109, 25, 53, 94,187,190,178,
         97, 49,126,  6, 64,143,169,173, 43,199, 33, 96, 76, 11,236, 23,
        166,104,141,246,125,217,122,238, 27, 22,228,146,165, 71,232,222
    };
    uint8_t h=seed;
    for(size_t i=0;i<nwords;i++) {
        h = Rand8[h ^ (uint8_t)(kmer[i]>>56)];
        h = Rand8[h ^ (uint8_t)(kmer[i]>>48) & 255];
        h = Rand8[h ^ (uint8_t)(kmer[i]>>40) & 255];
        h = Rand8[h ^ (uint8_t)(kmer[i]>>32) & 255];
        h = Rand8[h ^ (uint8_t)(kmer[i]>>24) & 255];
        h = Rand8[h ^ (uint8_t)(kmer[i]>>16) & 255];
        h = Rand8[h ^ (uint8_t)(kmer[i]>>8) & 255];
        h = Rand8[h ^ (uint8_t)kmer[i] & 255];
    }
    return h;
}

inline size_t kmerizer::pos2kmer(size_t pos, kword_t *kmer, vector<bvec*> &index) {
    if (pos >= index[0]->get_size()) return 1;
    size_t bpw = 8*sizeof(kword_t);
    // need to zero the kmer first
    memset(kmer,0,kmer_size);
    for(size_t b=0;b<8*kmer_size;b++)
        if (index[b]->find(pos))
            kmer[b/bpw] |= 1ULL << (bpw - b%bpw - 1);
    return 0;
}

inline uint32_t kmerizer::pos2value(size_t pos, vector<uint32_t> &values, vector<bvec*> &index) {
    // lookup the value in the pos bit
    // find the first bvec where this bit is set
    for(size_t i=0;i<values.size();i++)
        if (index[i]->find(pos))
            return values[i];
    return 0;
}


#endif // #ifndef _KMERIZER_H_
