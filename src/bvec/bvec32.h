#ifndef SNAPDRAGON_BVEC32_H
#define SNAPDRAGON_BVEC32_H

#include <stdint.h>

#define DEBUG false

#define WORD_SIZE               32
#define LITERAL_SIZE            31
#define BIT1          0x80000000UL
#define BIT2          0x40000000UL
#define FILLMASK      0x3FFFFFFFUL
#define ALL1S         0x7FFFFFFFUL
#define ONEFILL       0xC0000000UL
#define ONEFILL1      0xC0000001UL
#define ONEFULL       0xFFFFFFFFUL
#define ZEROFULL      0xAFFFFFFFUL
#define ZEROFILL1     0x80000001UL

typedef uint32_t word_t;

#endif // SNAPDRAGON_BVEC32_H