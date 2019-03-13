#ifndef cola_h
#define cola_h

#include <stdint.h>
#include <string.h>
#include <stdio.h>



void cola_keygen(
    int16_t *h,
    int16_t *f
);

void cola_encaps(
    const int16_t *h, 
    int16_t *c,
    unsigned char *k
);

void cola_decaps(
    int16_t *c,
    const int16_t *h,
    const int16_t *f,
    unsigned char *k
);
#endif