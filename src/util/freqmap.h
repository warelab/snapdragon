#ifndef SNAPDRAGON_FREQMAP_H
#define SNAPDRAGON_FREQMAP_H

#include <map>

class FrequencyMap {
    std::map<const char *, int> map;
public:
    void add(const char *key);
    int  count(const char *key);
};

#endif // #ifndef SNAPDRAGON_FREQMAP_H