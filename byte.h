#include <stdint.h>

//字节数组长度 = 735
void ZqArray2byteArray(const int64_t *aa, unsigned char *b);
void byteArray2ZqArray(const unsigned char *b, int64_t *aa);

/*588 个分量, 每个分量占2比特表示{-1,0,1}时, 每 4个系数, 占一个字节, 所以 588/4=147个字节
 */
//字节数组长度 = 147
void trinary2byteArray(const int64_t *f, unsigned char *b);
void byteArray2trinary(unsigned char *b, int64_t *f);

// 字节数组长度=74
void binary2byteArray(const int64_t *m, unsigned char *b);
void byteArray2binary(const unsigned char *b, int64_t *m);

void poly2bytes(unsigned char *b, const int64_t * f);

void bytes2poly(int64_t *f, const unsigned char *b);

void displaybits(unsigned char *b);

void display10bits(int64_t *f);