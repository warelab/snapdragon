#include "freqmap.h"

#include <map>

int FrequencyMap::count(const char *key) {
    return map[key];
}

void FrequencyMap::add(const char *key) {
    map[key]++;
}