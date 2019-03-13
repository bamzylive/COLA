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
    int16_t *h,
    int16_t *f
)
{
    printf("Z_1024[X] / X^%hd - 1\n",N_DEG);
    int16_t g[N], f_inv[N];
    memset( f_inv, 0, N*sizeof(int16_t) );
    memset( h, 0, N*sizeof(int16_t));
    srand(time(NULL));
    trinary_poly_gen(g,TRI_d);
    trinary_poly_gen(f,TRI_d);
    for(int i = 0; i < N; i++){
        g[i] += Q>>1 ;
    }
    int cnt = 20;
    while(!lift_power2_inv(f_inv,f,N_DEG) && cnt--){
        trinary_poly_gen(f,TRI_d);
    }
    cyc_convolution(h, f_inv, g, N_DEG);
    entrywise_mod_p(h, Q);
    
    central_mod_p(h, Q);
    // printf("pk is\n");
    if (dg(h)!=-1) printf("pk is OK.\n");
    else printf("pk is ZERO.\n");
    // display(h);
    // printf("sk is\n");
    // display(f);
}

void cola_encaps(
    const int16_t *h, 
    int16_t *c,
    unsigned char *k
)
{
    int i;
    int16_t r[N],e[N];
    binary_poly_gen(r);
    binary_poly_gen(e);
    cyc_convolution(c, h, r, N_DEG);
    for( i = 0; i<N; i++){
        c[i] += e[i];
    } 
    entrywise_mod_p(c,Q);
    central_mod_p(c, Q);
    // HASH(h,r,c)
    unsigned char str2hash[3*BYTESLEN];
    poly2bytes(str2hash, h);
    poly2bytes(str2hash+BYTESLEN, r);
    poly2bytes(str2hash+2*BYTESLEN, c);
    shake256(k, 256, str2hash, 3*BYTESLEN);
}

void cola_decaps(
    int16_t *c,
    const int16_t *h,
    const int16_t *f,
    unsigned char *k
)
{
    int i;
    int16_t d[N],r[N];
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
    /* int16_t h[N],f[N],c[N];
    unsigned char k1[256+1], k2[256+1];
    cola_keygen(h, f); 
    printf("still here?\n");
    cola_encaps(h, c, k1);
    cola_decaps(c, h, f, k2);
    k1[256] = 0;
    k2[256] = 0;
    if (strcmp(k1, k2) == 0) printf("succeeded\n");
    else printf("failed\n"); */
    
    int16_t a[N], b[N], s[N], t[N],d[N];
    memset(a, 0, sizeof(int16_t)*N);
    memset(b, 0, sizeof(int16_t)*N);
    int16_t p = 2;

    //a[4] = 1; a[3] = 3; a[1] = 2; a[0] = 4; 
    a[4] = a[2] = a[0] = 1;

    //b[2] = 1;  b[0] = -1;
    
    display(a);
    quotient_ring_inv(b,a,p,7);
    printf("Its inverse is\n");
    display(b);
}
