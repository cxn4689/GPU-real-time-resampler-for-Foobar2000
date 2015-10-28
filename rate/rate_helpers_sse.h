#if defined(SSE3_)
#define _mm_lddqu_ps(x) _mm_castsi128_ps(_mm_lddqu_si128((__m128i*)(x)))
#else
#define _mm_lddqu_ps(x) _mm_loadu_ps(x)
#endif


static __inline float _mm_add_horz(__m128 x)
{
#if defined(SSE3_)
    x = _mm_hadd_ps(x, x);
    x = _mm_hadd_ps(x, x);
#else
#pragma warning(disable : 4700)
    __m128 y;
    y = _mm_movehl_ps(y, x);
    x = _mm_add_ps(x, y);
#pragma warning(default : 4700)
    y = x;
    y = _mm_shuffle_ps(y, y, _MM_SHUFFLE(1,1,1,1));
    x = _mm_add_ss(x, y);
#endif
    return _mm_cvtss_f32(x);
}

static __inline __m128 _mm_add_horz_ss(__m128 x)
{
#if defined(SSE3_)
    x = _mm_hadd_ps(x, x);
    x = _mm_hadd_ps(x, x);
#else
    __m128 y;
#ifdef _DEBUG
    y = _mm_setzero_ps();
#endif
    y = _mm_movehl_ps(y, x);
    x = _mm_add_ps(x, y);
    y = x;
    y = _mm_shuffle_ps(y, y, _MM_SHUFFLE(1,1,1,1));
    x = _mm_add_ss(x, y);
#endif
    return x;
}


/* good idea to make them static __inline func() */


#define SUM1(n) \
        xmm0 = _mm_lddqu_ps(input+n);\
        xmm1 = _mm_lddqu_ps(input-n);\
        xmm2 = _mm_load_ps(coeff-n);\
        \
        xmm3 = _mm_move_ss(xmm3, xmm0);\
        xmm3 = _mm_shuffle_ps(xmm3, xmm3, _MM_SHUFFLE(1, 2, 3, 0));\
        xmm3 = _mm_add_ps(xmm3, xmm1);\
        xmm3 = _mm_mul_ps(xmm3, xmm2);\
        xmm7 = _mm_add_ps(xmm7, xmm3)


#define SUM2(n) \
        xmm3 = _mm_lddqu_ps(input+n);\
        xmm4 = _mm_lddqu_ps(input-n);\
        xmm5 = _mm_load_ps(coeff-n);\
        \
        xmm0 = _mm_move_ss(xmm0, xmm3);\
        xmm0 = _mm_shuffle_ps(xmm0, xmm0, _MM_SHUFFLE(1, 2, 3, 0));\
        xmm0 = _mm_add_ps(xmm0, xmm4);\
        xmm0 = _mm_mul_ps(xmm0, xmm5);\
        xmm7 = _mm_add_ps(xmm7, xmm0)



#define ADD0A() \
    xmm0 = _mm_load_ps(coef);\
    xmm2 = _mm_lddqu_ps(at);\
    xmm3 = _mm_load_ps(coef+4);\
    xmm4 = _mm_lddqu_ps(at+4);\
    xmm0 = _mm_mul_ps(xmm0, xmm2);\
    xmm3 = _mm_mul_ps(xmm3, xmm4);\
    xmm0 = _mm_add_ps(xmm0,xmm3)

#define ADD0B(n) \
    xmm1 = _mm_load_ps(coef+n);\
    xmm2 = _mm_lddqu_ps(at+n);\
    xmm3 = _mm_load_ps(coef+n+4);\
    xmm4 = _mm_lddqu_ps(at+n+4);\
    xmm1 = _mm_mul_ps(xmm1, xmm2);\
    xmm0 = _mm_add_ps(xmm0, xmm1);\
    xmm3 = _mm_mul_ps(xmm3, xmm4);\
    xmm0 = _mm_add_ps(xmm0,xmm3);\


#define ADD1p(n) \
    xmm1 = _mm_load_ps(coef+2*n);\
    xmm2 = _mm_load_ps(coef+2*n+4);\
    xmm4 = _mm_load_ps(coef+2*n+8);\
    xmm5 = _mm_load_ps(coef+2*n+12);\
    xmm3 = _mm_lddqu_ps(at+n);\
    xmm6 = _mm_lddqu_ps(at+n+4);\
    xmm1 = _mm_mul_ps(xmm1, x);\
    xmm2 = _mm_add_ps(xmm2, xmm1);\
    xmm3 = _mm_mul_ps(xmm3, xmm2);\
    xmm4 = _mm_mul_ps(xmm4, x);\
    xmm5 = _mm_add_ps(xmm5, xmm4);\
    xmm6 = _mm_mul_ps(xmm6, xmm5);\
    sum = _mm_add_ps(sum, xmm3);\
    sum = _mm_add_ps(sum, xmm6)



#define ADD2(n) \
    xmm1 = _mm_load_ps(coef+3*n);\
    xmm2 = _mm_load_ps(coef+3*n+4);\
    xmm3 = _mm_load_ps(coef+3*n+8);\
    xmm4 = _mm_lddqu_ps(at+n);\
    xmm1 = _mm_mul_ps(xmm1, x2);\
    xmm2 = _mm_mul_ps(xmm2, x);\
    xmm3 = _mm_add_ps(xmm3, xmm1);\
    xmm3 = _mm_add_ps(xmm3, xmm2);\
    xmm3 = _mm_mul_ps(xmm3, xmm4);\
    sum = _mm_add_ps(sum, xmm3)
