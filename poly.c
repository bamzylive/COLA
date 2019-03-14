#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include "fastrandombytes.h"
#include "params.h"
#include "poly.h"

void binary_poly_gen(int64_t *f)
{
    memset(f, 0, sizeof(int64_t) * N);
    uint64_t r;
    uint64_t i, j, index;
    for (i = 0; i <= N_DEG / 16; i++)
    {
        r = rand() & ((1 << 16) - 1);
        for (j = 0; j < 16; j++)
        {
            index = i * 16 + j;
            if (index < N_DEG)
                f[index] = (r & (1 << j)) >> j;
        }
    }
}

void trinary_poly_gen(
    int64_t *f,
    const uint64_t d)
{

    uint64_t r, i, count, coeff[6];

    memset(f, 0, sizeof(int64_t) * N); // string.h
    count = 0;
    while (count < d + 1)
    {
        rng_uint64(&r); // 64-bit random binary string
        for (i = 0; i < 6; i++)
        {
            coeff[i] = r & 0x3FF;     // least 10 bits remain the same, others 0
            r = (r - coeff[i]) >> 10; // chop off last 10 bits
            if (coeff[i] < N_DEG)
            {
                if (f[coeff[i]] == 0)
                {
                    f[coeff[i]] = 1;
                    count++;
                    if (count == d + 1)
                        break;
                }
            }
        }
    }
    count = 0;
    while (count < d)
    {
        rng_uint64(&r);
        for (i = 0; i < 6; i++)
        {
            coeff[i] = r & 0x3FF;
            r = (r - coeff[i]) >> 10;
            if (coeff[i] < N_DEG)
            {
                if (f[coeff[i]] == 0)
                {
                    f[coeff[i]] = -1;
                    count++;
                    if (count == d)
                        break;
                }
            }
        }
    }
    return;
}

void gen_rand_poly(int64_t *f, int64_t f_deg, int64_t num)
{
    for (int i = 0; i <= f_deg; i++)
        f[i] = (rand()) % num;
}

int64_t dg(const int64_t *f)
{
    int64_t degree = -1; // -1 represents zero polynomial
    for (int i = N - 1; i > -1; i--)
        if (f[i] != 0)
        {
            degree = i;
            break;
        }
    return degree;
}

int64_t display(const int64_t *f)
{
    int64_t i;
    int64_t degree = dg(f);
    if (degree < 0)
    {
        printf("[  0  ]  zero polynomial\n");
        return degree;
    }
    printf("[  ");
    for (i = 0; i <= degree; i++)
    {
        if (f[i])
            printf("( %lld * x^ %lld )+", f[i], i);
    }
    printf("] of degree %lld\n", degree);

    return degree;
}

/* non-zero poly assignment */
void poly_assign(int64_t *y, int64_t *x)
{

    for (int i = 0; i < N; i++)
    {
        y[i] = x[i];
    }
}

// output gcd
int64_t EEA(
    int64_t a,
    int64_t b,
    int64_t *ps,
    int64_t *pt)
{
    int64_t r0, r1, r2, s0, s1, s2, t0, t1, t2, q2;
    r0 = a;
    r1 = b;
    s0 = 1;
    s1 = 0;
    t0 = 0;
    t1 = 1;

    while (r0 % r1 != 0)
    {
        r2 = r0 % r1;
        q2 = (r0 - r2) / r1;
        s2 = s0 - q2 * s1;
        t2 = t0 - q2 * t1;
        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }
    *ps = s1;
    *pt = t1;
    return r1;
}

/* 
*   return a number inv in [0,p-1] 
*   p must be a prime  
*/
int64_t num_inv(int64_t lc, int64_t p)
{
    int64_t s, t;
    EEA(lc, p, &s, &t);
    return s < 0 ? s + p : s;
}

/*
  input: a,b, degree of a, degree of b, q and r, here q and r
  length of all poly is N
  output: q, r such that a = qb + r
*/
void poly_div(
    int64_t p,
    const int64_t *a,
    const int64_t *b,
    int64_t *q,
    int64_t *r,
    int64_t a_deg,
    int64_t b_deg)
{
    if (p == 2)
    {
        bin_poly_div(a, b, q, r, a_deg, b_deg);
        return;
    }

    int i, j;
    int64_t u = num_inv(b[b_deg], p); // inv of leading coef of b mod p
    int64_t r_deg, d;
    int64_t v[N];

    //   printf("The dividend is \n");
    //   display(a);
    //   printf("The divisor is \n");
    //   display(b);

    memset(q, 0, N * sizeof(int64_t));
    memset(r, 0, N * sizeof(int64_t));
    memset(v, 0, N * sizeof(int64_t));

    for (i = 0; i <= a_deg; i++)
    {
        r[i] = a[i];
    }

    r_deg = a_deg;
    //printf("initially degree of r is %lld\n", r_deg);
    //printf("initially degree of b is %lld\n", b_deg);
    while (r_deg >= b_deg)
    {
        //printf("current degree of r is %lld\n", r_deg);
        //printf("current degree of b is %lld\n", b_deg);
        d = r_deg;
        v[d - b_deg] = (u * r[d]) % p;
        for (j = r_deg; j >= 0; j--)
        {
            r[j] = (r[j] - v[d - b_deg] * b[j - d + b_deg]) % p;
            if (r[j] < 0)
                r[j] += p;

            if (r[j] != 0)
            {
                r_deg = j;
                //printf("now degree of r is DOWN to %hu\n", r_deg);
                break;
            }
        }
        if (j == -1)
        {
            r_deg = -1;
        }
        for (j = r_deg - 1; j >= 0; j--)
        {
            r[j] = (r[j] - b[j - d + b_deg] * v[d - b_deg]) % p;
            if (r[j] < 0)
                r[j] += p;
        }
        q[d - b_deg] = v[d - b_deg];
    }

    // printf("quotient is\n");
    // display(q);
    // printf("remainder is \n");
    // display(r);
}

void bin_poly_div(
    const int64_t *a,
    const int64_t *b,
    int64_t *q,
    int64_t *r,
    int64_t a_deg,
    int64_t b_deg)
{
    int i, j;
    int64_t r_deg, d;

    // printf("The dividend is \n");
    // display(a);
    // printf("The divisor is \n");
    // display(b);

    memset(q, 0, N * sizeof(int64_t));
    memset(r, 0, N * sizeof(int64_t));

    for (i = 0; i <= a_deg; i++)
    {
        r[i] = a[i] & 1; // 100000000 - |x| : -1 -> 11111111, -2 -> 11111110
    }

    r_deg = a_deg;
    //printf("initially degree of r is %lld\n", r_deg);
    //printf("initially degree of b is %lld\n", b_deg);
    while (r_deg >= b_deg)
    {
        //printf("current degree of r is %lld\n", r_deg);
        //printf("current degree of b is %lld\n", b_deg);
        d = r_deg;
        for (j = r_deg; j >= 0; j--)
        {
            r[j] ^= (b[j - d + b_deg] & 1);

            if (r[j] != 0)
            {
                r_deg = j;
                //printf("now degree of r is DOWN to %hu\n", r_deg);
                break;
            }
        }
        if (j == -1)
        {
            r_deg = -1;
        }
        for (j = r_deg - 1; j >= 0; j--)
        {
            r[j] ^= (b[j - d + b_deg] & 1);
        }
        q[d - b_deg] = 1;
    }
}

/*
  textbook mul over Z
*/
void textbook_mul(
    int64_t *c,
    int64_t *a,
    int64_t *b)
{
    int64_t b_deg = dg(b);
    int64_t a_deg = dg(a);
    int64_t d;
    memset(c, 0, N*sizeof(int64_t));

    if (b_deg * a_deg < 0)
        return;
    d = (a_deg > b_deg) ? a_deg : b_deg;
    int i, j;
    int64_t temp;
    for (i = 0; i <= d; i++)
    {
        temp = 0;
        for (j = 0; j <= i; j++)
            temp += a[j] * b[i - j];
        c[i] = temp;
    }
    for (i = d + 1; i <= a_deg + b_deg; i++)
    {
        temp = 0;
        for (j = d; j >= i - d; j--)
            temp += b[j] * a[i - j];
        c[i] = temp;
    }
}
/*
  karatsuba mul over Z
*/
void karatsuba(
    int64_t *c,
    int64_t *a,
    int64_t *b)
{
    int64_t b_deg = dg(b);
    int64_t a_deg = dg(a);
    int64_t d, n;
    printf("degree of a is %lld\n", a_deg);
    printf("degree of b is %lld\n", b_deg);
    if (b_deg * a_deg < 0)
    {
        printf("here\n");
        return;
    }
    d = (a_deg > b_deg) ? a_deg : b_deg;
    n = d + 1;

    int64_t dd[N];
    memset(dd, 0, N * sizeof(int64_t));
    for (int i = 0; i < n; i++)
    {
        dd[i] = a[i] * b[i];
    }
    c[0] = dd[0];
    if (n == 1)
        return;
    c[2 * n - 2] = dd[n - 1];
    int64_t sum1, sum2;
    for (int64_t i = 1; i < 2 * n - 2; i++)
    {
        sum1 = 0;
        sum2 = 0;
        for (int64_t p = 0; p <= (i - 1) / 2; p++)
        {
            sum1 += (dd[p] + dd[i - p]);
            sum2 += ((a[p] + a[i - p]) * (b[p] + b[i - p]));
        }
        if ((i & 1) == 1)
        {
            c[i] = sum2 - sum1;
        }
        else
        {
            c[i] = sum2 - sum1 + dd[i / 2];
        }
    }
}

void poly_minus(int64_t *c, int64_t *a, int64_t *b)
{
    for (int i = 0; i < N; i++)
    {
        c[i] = a[i] - b[i];
    }
}

void entrywise_mod_p(int64_t *c, int64_t p)
{
    // if (p == 0) return;
    if (p == 2)
    {
        for (int i = 0; i <= dg(c); i++)
            c[i] &= 1;
        return;
    }
    for (int i = 0; i <= dg(c); i++)
    {
        c[i] = c[i] % p;
        if (c[i] < 0)
            c[i] += p;
    }
}

void central_mod_p(int64_t *c, int64_t p)
{
    for (int i = 0; i <= dg(c); i++)
    {
        if (c[i] >= (p >> 1))
        {
            c[i] = c[i] - p;
        }
    }
}

void poly_EEA(
    int64_t *a,
    int64_t *b,
    int64_t *olds,  // size: 2*N
    int64_t *oldt,  // size: 2*N
    int64_t *oldr,  // size: 2*N
    int64_t p)
{
    int64_t a_deg, b_deg;
    a_deg = dg(a);
    b_deg = dg(b);
    memset(olds, 0, N*sizeof(int64_t));
    memset(oldt, 0, N*sizeof(int64_t));
    memset(oldr, 0, N*sizeof(int64_t));
    if (b_deg == -1) {
        olds[0] = 1;
        poly_assign(oldr, a);
        entrywise_mod_p(oldr, p);
        return;
    }
    if (a_deg == -1) {
        oldt[0] = 1;
        poly_assign(oldr, b);
        entrywise_mod_p(oldr, p);
        return;
    }
    int64_t s[N], t[N], r[N];
    int64_t quo[N], prov[N], rem[N], temp[N];
    memset(s, 0, N*sizeof(int64_t));
    olds[0] = 1;
    memset(t, 0, N*sizeof(int64_t));
    t[0] = 1;
    memset(r, 0, N*sizeof(int64_t));
    poly_assign(oldr, a);
    entrywise_mod_p(oldr, p);  // display(oldr);
    poly_assign(r, b);
    entrywise_mod_p(r, p);
    memset(quo, 0, N*sizeof(int64_t));

    while (dg(r) != -1) {
        memset(rem, 0, N*sizeof(int64_t));
        // quo = oldr / r
        poly_div(p, oldr, r, quo, rem, dg(oldr), dg(r));

        // (oldr, r) := (r, old_r - quotient*r)
        memset(prov, 0, N*sizeof(int64_t));
        poly_assign(prov, r);
        // printf("sum of degrees of quo and r is: %ld\n", dg(quo)+dg(r));
        textbook_mul(temp, quo, r);
        poly_minus(r, oldr, temp);
        entrywise_mod_p(r, p);
        poly_assign(oldr, prov);

        // (olds, s) := (s, old_s - quotient*s)
        memset(prov, 0, N*sizeof(int64_t));
        poly_assign(prov, s);
        // printf("sum of degrees of quo and s is: %ld\n", dg(quo)+dg(s));
        textbook_mul(temp, quo, s);
        poly_minus(s, olds, temp);
        entrywise_mod_p(s, p);
        poly_assign(olds, prov);

        // (oldt, t) := (t, old_t - quotient*t)
        memset(prov, 0, N*sizeof(int64_t));
        poly_assign(prov, t);
        // printf("sum of degrees of quo and t is: %ld\n", dg(quo)+dg(t));
        textbook_mul(temp, quo, t);
        poly_minus(t, oldt, temp);
        entrywise_mod_p(t, p);
        poly_assign(oldt, prov);
    }
    // output Bezout coef (olds, oldt)
    // gcd oldr
}

/*
 *    inverse in Zp[X] / ( X^n - 1 )
 *    N = n+1
 *    p must be a prime
 *    output: 1: a_inv s.t a * a_inv = 1 in Rp
 *            0: not invertible  
 */
int quotient_ring_inv(
    int64_t *a_inv,
    int64_t *a,
    int64_t p,
    int64_t n)
{
    // printf("Z_%hu[X] / X^%lld - 1\n",p,n);
    int64_t sum;
    int64_t a_deg = dg(a);
    memset(a_inv, 0, N * sizeof(int64_t));
    // printf("Whether a is a non-zero const\n");
    // a is a non-zero const
    if (a_deg == 0)
    {
        a_inv[0] = num_inv(a[0], p);
        return 1;
    }
    // printf("Whether monimial\n");
    // (monimial) if a = c_k X^k, then its INV is c_k^(-1) X^(n-k)
    sum = 0;
    for (int i = 0; i < a_deg; i++)
    {
        sum += a[i];
        if (sum != 0)
            break;
    }
    if (sum == 0)
    {
        a_inv[n - a_deg] = num_inv(a[a_deg], p);
        return 1;
    }
    //printf("mod 2 case: if sum of coefficients is even, NOT invertible\n");
    // mod 2 case: if sum of coefficients is even, NOT invertible
    if (p == 2)
    {
        sum = 0;
        for (int j = 0; j <= n + 1; j++)
            sum ^= (a[j] & 1);
        if (sum == 0)
            return 0;
    }
    // printf("poly_EEA routine\n");
    int64_t minus_one[N], u[N], v[N], gcd[N];
    memset(minus_one, 0, N * sizeof(int64_t));
    minus_one[0] = -(minus_one[n] = 1);
    // display(minus_one);
    poly_EEA(a, minus_one, u, v, gcd, p);
    // printf("u is\n");
    // display(u);
    // printf("v is\n");
    // display(v);
    // printf("gcd is\n");
    // display(gcd);
    // printf("");
    if (dg(gcd) == 0)
    {   
        // display(u);
        // int64_t prod[N];
        // memset(prod, 0, N*sizeof(int64_t));
        // cyc_convolution(prod, a, u, n);
        // entrywise_mod_p(prod, p);
        // display(prod);
        int64_t gcd_inv = num_inv(gcd[0], p);
        for (int i = 0; i < N; i++)
            a_inv[i] = gcd_inv * u[i] % p;
        //printf("here and now\n");
        return 1;
    }
    else
        return 0;
}

/*
    mul in Z[X] / ( X^n - 1 )
*/
void cyc_convolution(
    int64_t *c,
    const int64_t *a,
    const int64_t *b,
    int64_t n)
{
    memset(c, 0, N * sizeof(int64_t));
    int l, i, j;
    int64_t sum;

    for (l = 0; l < n; l++)
    {
        sum = 0;
        for (i = 1; i <= n; i++)
        {
            sum += b[(i + l) % n] * a[n - i];
        }
        c[l] = sum;
    }
}

/*
    To calculate INV in Z_2^r[X] / ( X^n - 1 ).  
*/
int lift_power2_inv(
    int64_t *b,
    int64_t *a,
    int64_t n,
    int64_t r)
{
    int i, j, flag;
    int64_t index = 1;
    int64_t temp1[N], temp2[N], modulus = 2;

    flag = quotient_ring_inv(b, a, 2, n);
    if (flag == 0)
    {
        printf("NOT INVERTIBLE!\n");
        return 0;
    }
    // display(a);
    // display(b);
    // int64_t prod[N];
    // memset(prod, 0, N*sizeof(int64_t));
    // cyc_convolution(prod, a, b, n);
    // entrywise_mod_p(prod, 2);
    // display(prod);

    while (index < r)
    {
        index *= 2;
        //printf("current index: %ld\n", index);
        modulus = modulus * modulus;
        //printf("current modulus: %ld\n", modulus);
        cyc_convolution(temp1, b, b, n);
        //display(temp1);
        entrywise_mod_p(temp1, modulus);
        cyc_convolution(temp2, temp1, a, n);
        entrywise_mod_p(temp2, modulus);
        for (i = 0; i < N; i++)
        {
            b[i] = b[i] * 2 - temp2[i];
        }
        entrywise_mod_p(b, modulus);
    }
    modulus = 1;
    for (i = 0; i < r; i++)
        modulus *= 2;
    //printf("final modulus: %ld\n", modulus);
    entrywise_mod_p(b, modulus);
    return 1;
}
