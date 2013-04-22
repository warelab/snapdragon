#ifndef SNAPDRAGON_BVEC_H
#define SNAPDRAGON_BVEC_H

#define DEBUG false

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cstring>
#include <boost/serialization/vector.hpp>
#include <boost/cstdint.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace std;

typedef uint32_t word_t;

class bvec {
    vector<word_t> words;
    bool rle;
    word_t count; // cache the number of set bits
    word_t size; // bits in the uncompressed bitvector

    friend class boost::serialization::access;
    friend ostream & operator <<(ostream &, const bvec &);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & words & count & size & rle;
    }

    // checkpoint: active word and associated info
    struct checkpoint {
        vector<word_t>::iterator active_word;
        word_t bit_pos;
    } frontier;

public:
    // Destructor
    ~bvec() {};

    // Constructors
    bvec() : rle(false), count(0), size(0) {};
    bvec(bool wah) : rle(false), count(0), size(0) {compress();};
    bvec(vector<word_t>& vals);
    bvec(word_t* buf); // DIY deserialization
    void print();
    void compress();
    void decompress();
    size_t dump(word_t **buf); // DIY serialization

    // logical set operations
    void flip();
    bvec* copyflip();
    void operator|=(bvec& rhs);
    bvec* operator|(bvec&);
    void operator&=(bvec& rhs);
    bvec* operator&(bvec&);
    bool operator==(bvec&) const;
    
    bool equals(const bvec&) const;
    vector<word_t>& get_words();

    // is x in the set?
    bool find(word_t x);
    // find the position of the next set bit after x
    word_t next_one(word_t x);
    
    // insert x into an existing bvec (at the end is faster)
    void setBit(word_t x);
    // for constructing a rle bvec one bit at a time
    void appendFill(bool bit, word_t count);

    // basic metrics
    word_t cnt();
    word_t get_size();
    word_t bytes();
    
    static void save_to_file(const bvec &bv, const char *filename);
    static void restore_from_file(bvec &bv, const char *filename);

private:

    bool low_density(vector<word_t>& vals);
    void construct_rle(vector<word_t>& vals);
    bvec& copy(const bvec& bv);
    void matchSize(bvec& bv);
    void rle_OR_rle(bvec& rhs);
    void rle_OR_non(bvec& rhs);
    void non_OR_rle(bvec& rhs);
    void non_OR_non(bvec& rhs);
    void rle_AND_rle(bvec& rhs);
    void rle_AND_non(bvec& rhs);
    void non_AND_rle(bvec& rhs);
    void non_AND_non(bvec& rhs);
    word_t popcount(word_t val) const;

};

#endif
