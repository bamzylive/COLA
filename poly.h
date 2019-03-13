#ifndef poly_h
#define poly_h





void
binary_poly_gen(int16_t  *f);

/*
  T(d+1,d)
*/
void
trinary_poly_gen(
          int16_t  *f,
    const uint16_t  d);

void gen_rand_poly(int16_t *f,int16_t f_deg, int16_t num);

int16_t dg(const int16_t  *f);

int16_t display(const int16_t  *f);

void poly_assign(int16_t *y, int16_t *x);

int16_t  EEA(
  int16_t a,
  int16_t b,
  int16_t* ps,
  int16_t* pt);

int16_t num_inv(int16_t lc, int16_t p);

/*
  input: a,b, degree of a, degree of b, q and r, here q and r
  length of all poly is N
  output: q, r such that a = qb + r
*/
void poly_div(
    int16_t p,
   const int16_t  *a,
   const int16_t  *b,
   int16_t  *q,
   int16_t  *r,
   int16_t a_deg,
   int16_t b_deg
);

void bin_poly_div(
   const int16_t  *a,
   const int16_t  *b,
   int16_t  *q,
   int16_t  *r,
   int16_t a_deg,
   int16_t b_deg
);

/*
  textbook mul over Z
*/
void textbook_mul(
    int16_t *a,
    int16_t *b,
    int16_t *c
);

/*
  karatsuba mul over Z
*/
void karatsuba(
    int16_t *a,
    int16_t *b,
    int16_t *c
);

void poly_minus(int16_t* c,int16_t *a, int16_t*b);
void entrywise_mod_p(int16_t* c, int16_t p);
void central_mod_p(int16_t* c, int16_t p);

void poly_EEA(
  int16_t  *a,
  int16_t  *b,
  int16_t  *olds,
  int16_t  *oldt,
  int16_t  *oldr,
  int16_t p
);

/*
 *    inverse in Zp[X]/ ( X^n - 1 )
 *    N = n+1
 *    p must be a prime
 *    output: 1: a_inv s.t a * a_inv = 1 in Rp
 *            0: not invertible  
 */
int quotient_ring_inv(int16_t * a_inv, int16_t *a, int16_t p, int16_t n);

void cyc_convolution(
    int16_t *c, 
    const int16_t *a, 
    const int16_t *b,
    int16_t n);

int lift_power2_inv(
    int16_t *b,
    int16_t *a,
    int16_t n,
    int16_t r
);

#endif /* poly_h */
