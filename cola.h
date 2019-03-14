#ifndef cola_h
#define cola_h

#include <stdint.h>
#include <string.h>
#include <stdio.h>



void cola_keygen(
    int64_t *h,
    int64_t *f
);

void cola_encaps(
    const int64_t *h, 
    int64_t *c,
    int64_t *r,
    unsigned char *k
);

void cola_decaps(
    int64_t *c,
    int64_t *r,
    const int64_t *h,
    const int64_t *f,
    unsigned char *k
);
#endif