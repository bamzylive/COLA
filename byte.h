#include <stdint.h>

void poly2bytes(unsigned char *b, const int16_t * f);

void bytes2poly(int16_t *f, const unsigned char *b);

void displaybits(unsigned char *b);

void display10bits(int16_t *f);