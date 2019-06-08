#ifndef poly_h
#define poly_h




/* 生成系数在0,1里面随机取的多项式, 为了后面乘法运算类型一致, 用int64_t来存储系数 */
void
binary_poly_gen(int64_t  *f);

/*
  生成三值多项式
  T(d+1,d)
*/
void
trinary_poly_gen(
          int64_t  *f,
    const uint64_t  d);

void gen_rand_poly(int64_t *f,int64_t f_deg, int64_t num);

/**
 * 返回多项式的次数
 */
int64_t dg(const int64_t  *f);

int64_t display(const int64_t  *f);

void poly_assign(int64_t *y, int64_t *x);

/**
 * 多项式扩展的欧几里得算法
 */
int64_t  EEA(
  int64_t a,
  int64_t b,
  int64_t* ps,
  int64_t* pt);

int64_t num_inv(int64_t lc, int64_t p);

/*
  input: a,b, degree of a, degree of b, q and r, here q and r
  length of all poly is N
  output: q, r such that a = qb + r
*/
void poly_div(
    int64_t p,
   const int64_t  *a,
   const int64_t  *b,
   int64_t  *q,
   int64_t  *r,
   int64_t a_deg,
   int64_t b_deg
);

void bin_poly_div(
   const int64_t  *a,
   const int64_t  *b,
   int64_t  *q,
   int64_t  *r,
   int64_t a_deg,
   int64_t b_deg
);

/*
  textbook mul over Z
*/
void textbook_mul(
    int64_t *a,
    int64_t *b,
    int64_t *c
);

/*
  karatsuba mul over Z
*/
void karatsuba(
    int64_t *a,
    int64_t *b,
    int64_t *c
);

void poly_minus(int64_t* c,int64_t *a, int64_t*b);
void entrywise_mod_p(int64_t* c, int64_t p);
void central_mod_p(int64_t* c, int64_t p);

void poly_EEA(
  int64_t  *a,
  int64_t  *b,
  int64_t  *olds,
  int64_t  *oldt,
  int64_t  *oldr,
  int64_t p
);

/*
 *    inverse in Zp[X]/ ( X^n - 1 )
 *    N = n+1
 *    p must be a prime
 *    output: 1: a_inv s.t a * a_inv = 1 in Rp
 *            0: not invertible  
 */
int quotient_ring_inv(int64_t * a_inv, int64_t *a, int64_t p, int64_t n);

void cyc_convolution(
    int64_t *c, 
    const int64_t *a, 
    const int64_t *b,
    int64_t n);

int lift_power2_inv(
    int64_t *b,
    int64_t *a,
    int64_t n,
    int64_t r
);

#endif /* poly_h */
