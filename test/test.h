#ifndef SNAPDRAGON_TEST_H
#define SNAPDRAGON_TEST_H

#include <vector>
#include <stdio.h>
#include <stdint.h>
#include "gtest/gtest.h"

using namespace std;

#define TERM_RESET "\e[m"
#define TERM_BOLD  "\e[1m"
#define TERM_RED   "\e[31m"
#define TERM_GREEN "\e[32m"

void print_binary(vector<uint32_t>& v) {
    for (int i = 0; i < v.size(); i++) {
        for (int b = 31; b >= 0; b--)
            printf("%d", v[i] & (1 << b) ? 1 : 0);
        printf(" ");
    }
    printf("\n");
}

void xor_vectors(const vector<uint32_t>& a, const vector<uint32_t>& b,
    vector<uint32_t>& c) {
    for (int i = 0; i < a.size(); i++) {
        c.push_back(a[i] ^ b[i]);
    }
}

void debug_binary(const char* head, vector<uint32_t>& v) {
    printf("%-20s ", head);
    print_binary(v);
}

void print_test(const char * name, const bool passed) {
    printf(""TERM_BOLD"Test %-20s"TERM_RESET": %s\n", name, passed
        ? ""TERM_GREEN"Pass"TERM_RESET"" : ""TERM_RED"Fail"TERM_RESET"");    
}

#endif // #ifndef SNAPDRAGON_TEST_H