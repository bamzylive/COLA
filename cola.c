#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "cola.h"
#include "poly.h"
#include "byte.h"
#include "params.h"
#include "fips202.h"
#include "api.h"

uint64_t cola_keygen(
    int64_t *h,
    int64_t *f
)
{
    
    uint64_t nounce;
    int64_t g[N], f_inv[N];
    memset( f_inv, 0, N*sizeof(int64_t) );
    memset( h, 0, N*sizeof(int64_t));

    srand(time(NULL));

    trinary_poly_gen(g,TRI_d);
    nounce=trinary_poly_gen(f,TRI_d);
    
    g[0] += Q / 2;
    int cnt = 20;
    while(!lift_power2_inv(f_inv,f,N_DEG,10) && cnt--){
        trinary_poly_gen(f,TRI_d);
    }
    cyc_convolution(h, f_inv, g, N_DEG);
    entrywise_mod_p(h, Q);
    printf("The degree of public key h is of %lld\n", dg(h));
   /* if (dg(h)!=-1) printf("pk is OK.\n");
    else printf("pk is ZERO.\n"); */
    return nounce;
}

void cola_encaps(
    const int64_t *h, 
    int64_t *c,
    int64_t *r,
    unsigned char *k
)
{
    int i;
    int64_t e[N];
    binary_poly_gen(r);
    binary_poly_gen(e);
    cyc_convolution(c, h, r, N_DEG);
    for( i = 0; i<N; i++){
        c[i] += e[i];
    } 
    entrywise_mod_p(c,Q);

    // central_mod_p(c, Q);
    // HASH(h,r,c)
    unsigned char str2hash[3*BYTESLEN];
    poly2bytes(str2hash, h);
    poly2bytes(str2hash+BYTESLEN, r);
    poly2bytes(str2hash+2*BYTESLEN, c);
    shake256(k, 256, str2hash, 3*BYTESLEN);
}
/* 
int pke_enc(
    unsigned char *pk,
    unsigned char *m, 
    unsigned long long mlen,
    unsigned char *c, 
    unsigned long long *clen)
{
}
 */
void cola_enc(
    const int64_t *h, 
    int64_t *c,
    int64_t *c2,
    int64_t *r, // 用来对比恢复前和恢复后 秘密多项式 是否一致
    int64_t *m)
{
    int i;
    int64_t e[N];
    binary_poly_gen(r);
    binary_poly_gen(e);
    cyc_convolution(c, h, r, N_DEG);
    for( i = 0; i<N; i++){
        c[i] += e[i];
        c2[i] = m[i]^r[i];  //c2 生成完毕
    } 
    entrywise_mod_p(c,Q); //c 生成完毕
    
}

void cola_dec(
    const int64_t *c,
    const int64_t *c2,
    const int64_t *f,
    int64_t *r,
    int64_t *m
)
{
    int i;
    int64_t d[N];
    memset(d, 0, sizeof(d));
    memset(r, 0, sizeof(int64_t)*N);
    cyc_convolution(d,f,c,N_DEG);
    entrywise_mod_p(d, Q);
    central_mod_p(d, Q);
    for(i=0; i<N; i++){
        if (d[i]>= - Q/4 && d[i]< Q/4)  r[i] = 0;
        else r[i] = 1;
        m[i] = r[i]^c2[i];
    }


}


void cola_decaps(
    int64_t *c,
    int64_t *r,
    const int64_t *h,
    const int64_t *f,
    unsigned char *k
)
{
    int i;
    int64_t d[N];
    memset(d, 0, sizeof(d));
    memset(r, 0, sizeof(int64_t)*N);

    cyc_convolution(d,f,c,N_DEG);
    entrywise_mod_p(d, Q);
    central_mod_p(d, Q);
    for(i=0; i<N; i++){
        if (d[i]>= - Q/4 && d[i]< Q/4)  r[i] = 0;
        else r[i] = 1;
    }
    // TODO: HASH(h,r,c)
    unsigned char str2hash[3*BYTESLEN];
    poly2bytes(str2hash, h);
    poly2bytes(str2hash+BYTESLEN, r);
    poly2bytes(str2hash+2*BYTESLEN, c);
    shake256(k, 256, str2hash, 3*BYTESLEN);
}


int pke_keygen(
    unsigned char *pk, 
    unsigned char *sk)
{
    int64_t f[N], h[N];
    memset(f,0,sizeof(f));
    memset(h,0,sizeof(h));
    cola_keygen(h,f);
    return 0;
}

int pke_enc(
    unsigned char *pk,
    unsigned char *m, 
    unsigned long long mlen, //74
    unsigned char *c, 
    unsigned long long clen)   // clen = 735 + 74
{
    int64_t h[N],message[N], r[N], c1[N], c2[N];
    unsigned char byte_c1[735];
    unsigned char byte_c2[74];

    byteArray2ZqArray(pk,h);
    byteArray2binary(m,message);
    //735+74

    cola_enc(h,c1,c2,r,message);
    ZqArray2byteArray(c1,byte_c1);
    binary2byteArray(c2,byte_c2);   

    unsigned long long i; 

    for(i=0; i<735; i++)
    {
        c[i] = byte_c1[i];
    }
    for(i=735; i<clen; i++)
    {
        c[i-735] = byte_c2[i];
    }
    return 0;
}

int pke_dec(
    unsigned char *sk,
    unsigned char *c, 
    unsigned long long clen,
    unsigned char *m, 
    unsigned long long mlen)
{
    int64_t f[N],c1[N],c2[N],message[N],r[N];
    byteArray2trinary(sk,f);
    unsigned char byte_c1[735];
    unsigned char byte_c2[74];
    int i;
    for(i=0; i<735; i++)
    {
        byte_c1[i] = c[i];
    }
    for(i=0; i<74; i++)
    {
        byte_c2[i] = c[735+i];
    }
    byteArray2ZqArray(byte_c1,c1);
    byteArray2binary(byte_c2,c2);
    cola_dec(c1,c2,f,r,message);
    binary2byteArray(message,m);

    return 0;

}



int main(){
    FILE *fpt;
    fpt = fopen("/Users/xty/Desktop/COLA/info.dat","wb+");
    int64_t nounce[N];

    uint32_t arr[7]={735,147,74,74,74,74,74}; 
    //公钥长度, 私钥长度, 明文1长度, 明文2长度, 密文1长度, 密文2长度, 随机数长度
    fwrite(arr, sizeof(arr), 1, fpt);

    int64_t h[N], f[N], c[N], r1[N], r2[N], diff[N];
    int64_t m[N], m1[N], s1[N], s2[N], c2[N]; 
    unsigned char byte_enc_m[74];
    unsigned char byte_dec_m[74];
    unsigned char byte_enc_c[735];
    unsigned char byte_enc_c2[74];
    unsigned char byte_nounce[74];
    //unsigned char byte_dec_c[735];
    unsigned char b[735];
    unsigned char bb[147];

    
    int64_t h2[N],f2[N];

    int i,flag;
    flag = 1;
    

    unsigned char k1[256+1], k2[256+1];

    clock_t start,finish;
    double total_time;
    int loop = 1;
     
    start=clock();
    printf("所选的多项式环为: Z_1024[X] / X^%d - 1\n",N_DEG);
    int t;
    for(t=0; t<loop; t++){

    
        cola_keygen(h, f); 
        //接下来写入公钥和私钥, 先转成字节数组
        ZqArray2byteArray(h,b);
        fwrite(b,sizeof(b), 1, fpt); //写公钥字节数组
        trinary2byteArray(f,bb);
        fwrite(bb,sizeof(bb),1,fpt); //写私钥字节数组



        binary_poly_gen(m);
        cola_enc(h,c,c2,s1,m); //加密明文m
        cola_dec(c,c2,f,s2,m1); //解密得到明文m1
        poly_compare(m,m1);
        return 0;

       
        binary2byteArray(m,byte_enc_m);
        fwrite(byte_enc_m,sizeof(byte_enc_m),1,fpt); //写明文1(加密测试)字节数组
        fwrite(byte_dec_m,sizeof(byte_dec_m),1,fpt); //写明文2(解密测试)字节数组

        ZqArray2byteArray(c,byte_enc_c);
        binary2byteArray(c2,byte_enc_c2);
        fwrite(byte_enc_c,sizeof(byte_enc_c),1,fpt);
        fwrite(byte_enc_c2,sizeof(byte_enc_c2),1,fpt); //写密文1(解密测试)字节数组
        fwrite(byte_enc_c,sizeof(byte_enc_c),1,fpt);
        fwrite(byte_enc_c2,sizeof(byte_enc_c2),1,fpt); //写密文2(解密测试)字节数组

                                                        //这两者不是一致的吗?
                                                        
        binary2byteArray(s1,byte_nounce);
        fwrite(byte_nounce,sizeof(byte_nounce),1,fpt);
        fclose(fpt);

        

        /*
        cola_encaps(h, c, r1, k1);
        cola_decaps(c, r2, h, f, k2);
        memset(diff, 0, N*sizeof(int64_t));
        poly_minus(diff, r1, r2);
        int64_t cnt = 0;
        for (int64_t i = 0; i < N; i++) {
            if (diff[i]) {
                printf("%lld: %lld\n", i, diff[i]);
                cnt++;
            }
        }
        printf("%lld\n", cnt); */
    }
    finish=clock();
    total_time=(double)(finish-start)/CLOCKS_PER_SEC;
    printf("Totoal time %f ms\n", (total_time*1000/loop) );



}
