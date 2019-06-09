#ifndef cola_h
#define cola_h

#include <stdint.h>
#include <string.h>
#include <stdio.h>



uint64_t cola_keygen(
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

void cola_enc(
    const int64_t *h, 
    int64_t *c,
    int64_t *c2,
    int64_t *r, // 用来对比恢复前和恢复后 秘密多项式 是否一致
    int64_t *m);

void cola_dec(
    const int64_t *c,
    const int64_t *c2,
    const int64_t *h,
    const int64_t *f,
    int64_t *s,
    int64_t *m
);
#endif