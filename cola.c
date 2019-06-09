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
    int64_t *m
)
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
    const int64_t *h,
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

int main(){
    FILE *fpt;
    fpt = fopen("/Users/xty/Desktop/COLA/info.dat","wb+");
    uint64_t nounce;

    uint32_t arr[7]={735,147,80,80,80,80,8}; 
    //公钥长度, 私钥长度, 明文1长度, 明文2长度, 密文1长度, 密文2长度, 随机数长度
    fwrite(arr, sizeof(arr), 1, fpt);

    int64_t h[N], f[N], c[N], r1[N], r2[N], diff[N];
    int64_t m[N], m1[N], s1[N], s2[N], c2[N]; 
    unsigned char byte_m[74];
    unsigned char byte_m1[74];
    unsigned char b[735];
    unsigned char bb[147];

    memset(b,0,sizeof(b));
    memset(bb,0,sizeof(bb));
    memset(byte_m,0,sizeof(byte_m));

    int64_t h2[N],f2[N];

    int i,flag;
    flag = 1;
    binary_poly_gen(m);
    binary2byteArray(m,byte_m);
    binary_poly_gen(m1);
    poly_compare(m,m1);

    //byteArray2binary(byte_m,m1);
    //poly_compare(m,m1);

    return 0;

    unsigned char k1[256+1], k2[256+1];

    clock_t start,finish;
    double total_time;
    int loop = 1;
     
    start=clock();
    printf("所选的多项式环为: Z_1024[X] / X^%d - 1\n",N_DEG);
    for(int t=0; t<loop; t++){

    
        nounce = cola_keygen(h, f); 
        //接下来写入公钥和私钥, 先转成字节数组
        ZqArray2byteArray(h,b);
        fwrite(b,sizeof(b), 1, fpt);
        trinary2byteArray(f,bb);
        fwrite(bb,sizeof(bb),1,fpt);

       /*  
        byteArray2trinary(bb,f2);
        poly_compare(f,f2);
        */ 
       /* 
        byteArray2ZqArray(b,h2);
        poly_compare(h,h2); */



        /* 
        cola_enc(h,c,c2,s1,m);
        cola_dec(c,c2,h,f,s2,m);
        for(i=0; i<N; i++){
            if(s1[i]!=s2[i]){
                printf("failed decryption at position %d.\n",i);
                flag = 0;
                break;
            }
        }
        if(flag) printf("Successful decryption!\n");
        */

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
