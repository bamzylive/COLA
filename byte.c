#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include "params.h"

// 588 个分量(最高分量为0), 需要735个字节, 每5个字节 确定 4个系数
void ZqArray2byteArray(const int64_t *aa, unsigned char *b)
{
    u_int16_t a[N];
    for(int k=0; k<N;k++){
        if(aa[k]<0) a[k] = Q + aa[k];
        else a[k] = aa[k];
    }
    int i,kuai;
    u_int16_t duo=0;
    //要存734个字节
    for(i = 0; i < N ; i = i+4 ) // 每 4个系数 存成 5个字节
    {

        kuai = i/4;
        // 5*kuai + j
        b[5*kuai] = a[i]&255; //255: 取低8比特

        duo = (a[i]>>8)&3; //8:截断低8比特, 3:取高2比特
        b[5*kuai+1] = duo + ( (a[i+1]&63)<<2 ); //63: 取下一个块的低6比特, 左移2比特
        
        duo = (a[i+1]>>6)&15; //6: 截断低6比特, 15: 取高4比特
        b[5*kuai+2] = duo + ( (a[i+2]&15)<<4 ); //15: 取下一个块的低4比特, 左移4比特,

        duo = (a[i+2]>>4)&63; //4: 截断低4比特, 63: 取高6比特
        b[5*kuai+3] = duo + ( (a[i+3]&3)<<6 ); //3: 取下一个块的低2比特, 左移6比特

        duo = (a[i+3]>>2)&255; //255: 截断低2比特
        b[5*kuai+4] = duo; 

    }

}

void byteArray2ZqArray(const unsigned char *b, int64_t *aa)
{
    
    u_int16_t a[N];

    memset(a,0,sizeof(a));

    int i,kuai;
    u_int16_t duo=0;
    // 每 5个字节 还原 4个系数 734
    for(i=0; i<735; i=i+5)
    {
        kuai = i/5;
        //4*kuai + j
        a[4*kuai] = b[i]+ ( (b[i+1]&3)<<8 ); //1byte+ next byte low 2 bit, left 6 bit

        duo = b[i+1]>>2; // take 6 bit

        a[4*kuai+1] = duo + ( (b[i+2]&15)<<6 );//next byte low 4 bit, left 4 bit

        duo = b[i+2]>>4; // take 4 bit
        a[4*kuai+2] = duo + ( (b[i+3]&63)<<4 ); //next low 6 bit, left 2 bit

        duo = b[i+3]>>6; // take 2 bit
        a[4*kuai+3] = duo + ( b[i+4]<<2 );
    }
   
    for(i=0; i<N; i++) aa[i] = a[i];
}



//588 个分量, 每个分量占2比特表示{-1,0,1}时, 每 4个系数, 占一个字节, 所以 588/4=147个字节
void trinary2byteArray(const int64_t *f, unsigned char *b)
{
    int i, kuai;

    uint16_t a[N];
    memset(a,0,sizeof(a));
    for(i=0; i<N; i++){
        if(f[i] < 0)  a[i]= 3+f[i]; 
        else a[i] = f[i]; //-1变为2
    }
    
    for(i=0; i<N; i=i+4)
    {
        kuai = i/4;
        b[kuai] = a[i];
        b[kuai] += ((a[i+1]&3)<<2);
        b[kuai] += ((a[i+2])<<4);
        b[kuai] += ((a[i+3]<<6));
    }

    
}
void byteArray2trinary(unsigned char *b, int64_t *f)
{
    memset(f,0,sizeof(int64_t)*N);
    unsigned char temp;
    int i;
    for(i=0; i<147; i++)
    {
        temp = b[i];
        for(int j=0; j<4; j++){
            f[4*i+j] = temp & 3;
            if (f[4*i+j] == 2) f[4*i+j]=-1;
            temp = (temp>>2);
        }
    }
}

void binary2byteArray(const int64_t *m, unsigned char *b)
{
    /* for(int i=0; i<8; i++)
    {
        printf("%lld\t",m[i]);
    } */
    //printf("\n");
    uint16_t temp;
    int i;
    for(i=0; i<73; i++)
    {   temp=0;
        for(int j=0; j<8; j++)
        {
            b[i] += (m[8*i+j]<<temp);
            temp++; 
        }
    }
    //printf("%u\n",b[0]);
    b[73] = m[586]*4+m[585]*2+m[584];

}

void byteArray2binary(const unsigned char *b, int64_t *m)
{
    memset(m,0,sizeof(int64_t)*N);
    m[586] = (b[73]>>2);
    m[585] = (b[73]>>1)&1;
    m[584] = (b[73]&1);
    unsigned char temp;
    int i;
    for(i=0; i<73; i++)
    {
        temp = b[i];
        for(int j=0; j<8; j++)
        {
            m[8*i+j] = temp & 1;
            temp =( temp>>1);
        }
    }
    
}

void poly2bytes(unsigned char *b, const int64_t * f) {
    memset(b, 0, BYTESLEN);
    int64_t spacebias = 0, spaceleft = 8, inputleft, copylen;
    int i,j;
    for (i = 0, j = 0; i < N; i++) {
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
    int i,j;
    for (i = 0, j = 0; i < BYTESLEN; i++) {
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
    int i;
    for ( i = 0; i < BYTESLEN; i++ ) {
        for (int j = 0; j < 8; j++)
            printf("%d", (b[i]&(1<<j)) >> j);
        printf("  ");
    }
    printf("\n");
}

void display10bits(int64_t *f) {
    int i;
    for ( i = 0; i < N; i++ ) {
        for (int j = 0; j < 10; j++)
            printf("%lld", (f[i]&(1<<j)) >> j);
        printf("  ");
    }
    printf("\n");
}
