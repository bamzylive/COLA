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

void cola_keygen(
    int64_t *h,
    int64_t *f
)
{
    printf("Z_1024[X] / X^%d - 1\n",N_DEG);
    int64_t g[N], f_inv[N];
    memset( f_inv, 0, N*sizeof(int64_t) );
    memset( h, 0, N*sizeof(int64_t));
    srand(time(NULL));
    trinary_poly_gen(g,TRI_d);
    trinary_poly_gen(f,TRI_d);
    // for(int i = 0; i < N; i++){
    //     g[i] += Q>>1 ;
    // }
    g[0] += Q / 2;
    int cnt = 20;
    while(!lift_power2_inv(f_inv,f,N_DEG,10) && cnt--){
        trinary_poly_gen(f,TRI_d);
    }
    cyc_convolution(h, f_inv, g, N_DEG);
    entrywise_mod_p(h, Q);
    // central_mod_p(h, Q);
    // printf("pk is\n");
    if (dg(h)!=-1) printf("pk is OK.\n");
    else printf("pk is ZERO.\n");

    // int64_t prod[N];
    // memset(prod, 0, N*sizeof(int64_t));
    // cyc_convolution(prod, f, f_inv, N_DEG);
    // entrywise_mod_p(prod, Q);
    // display(prod);

    // display(h);
    // printf("sk is\n");
    // display(f);
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
    memset(d, 0, sizeof(N));
    memset(r, 0, sizeof(N));

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
    int64_t h[N], f[N], c[N], r1[N], r2[N], diff[N];
    unsigned char k1[256+1], k2[256+1];

    clock_t start,finish;
    double total_time;
    int loop = 10;
     
    start=clock();

    for(int t=0; t<loop; t++){

    
        cola_keygen(h, f); 
        //printf("still here?\n");
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
        printf("%lld\n", cnt);
    }
    finish=clock();
    total_time=(double)(finish-start)/CLOCKS_PER_SEC;
    printf("Totoal time %f ms\n", (total_time*1000/loop) );



}
