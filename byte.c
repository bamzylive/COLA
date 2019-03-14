#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "params.h"

void poly2bytes(unsigned char *b, const int64_t * f) {
    memset(b, 0, BYTESLEN);
    int64_t spacebias = 0, spaceleft = 8, inputleft, copylen;
    for (int i = 0, j = 0; i < N; i++) {
        //printf("copying the no.%d integer\n", i);
        inputleft = 10;
        while (inputleft) {
            copylen = spaceleft < inputleft ? spaceleft : inputleft;
            b[j] |= (f[i] & (((1<<copylen)-1)<<(10-inputleft))) >> (10-inputleft) << spacebias;
            inputleft -= copylen;
            //printf("%d bits copied, %d left\n", copylen, inputleft);
            spacebias += copylen;
            if (spacebias >= 8) {
                j++;
                //printf("%d rooms used\n", j);
                spacebias -= 8;
            }
            spaceleft = 8 - spacebias;
        }
    }
}

void printbits(int64_t l, int64_t c) {
    for (int j = 0; j < l; j++)
        printf("%lld", (c&(1<<j)) >> j);
    printf("  ");
}

void bytes2poly(int64_t *f, const unsigned char *b) {
    memset(f, 0, N*sizeof(int64_t));
    int64_t spacebias = 0, spaceleft = 10, inputleft, copylen;
    for (int i = 0, j = 0; i < BYTESLEN; i++) {
        //printf("\ncopying the no.%d byte\n", i);
        inputleft = 8;
        while (inputleft) {
            //printf("from no.%d bit ", 8-inputleft);
            //printf("to room %d, start at %d\n", j, spacebias);
            copylen = spaceleft < inputleft ? spaceleft : inputleft;
            f[j] |= (b[i] & (((1<<copylen)-1)<<(8-inputleft))) >> (8-inputleft) << spacebias;
            inputleft -= copylen;
            //printf("copied: "); printbits(copylen, b[i] & (((1<<copylen)-1)<<(8-inputleft)));
            //printf("%d bits copied, %d left\n", copylen, inputleft);
            spacebias += copylen;
            if (spacebias >= 10) {
                j++;
                //printf("%d rooms used\n", j);
                //printbits(10, f[j-1]);
                spacebias -= 10;
            }
            spaceleft = 10 - spacebias;
        }
    }
}

void displaybits(unsigned char *b) {
    for (int i = 0; i < BYTESLEN; i++) {
        for (int j = 0; j < 8; j++)
            printf("%d", (b[i]&(1<<j)) >> j);
        printf("  ");
    }
    printf("\n");
}

void display10bits(int64_t *f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 10; j++)
            printf("%lld", (f[i]&(1<<j)) >> j);
        printf("  ");
    }
    printf("\n");
}
