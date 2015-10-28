/* good idea to make them static __inline func() */

#define _mm_lddqu_pd(x) _mm_castsi128_pd(_mm_lddqu_si128((__m128i*)(x)))

#define SUM02(n, _mm_loadfun) \
			XMM1 = _mm_load_pd(coeff-n);\
			XMM2 = _mm_loadfun(input-n);\
			XMM0 = _mm_load_sd(input+n);\
			XMM3 = _mm_loadfun(input+n-2); /*out XMM3*/\
			XMM0 = _mm_shuffle_pd(XMM0, XMM3, _MM_SHUFFLE2(1, 0));\
			XMM0 = _mm_add_pd(XMM0, XMM2);\
			XMM0 = _mm_mul_pd(XMM0, XMM1)

#define SUM1(n, _mm_loadfun) \
			XMM1 = _mm_load_pd(coeff-n);\
			XMM5 = _mm_load_pd(coeff-n+2);\
			XMM2 = _mm_loadfun(input-n);\
			XMM6 = _mm_loadfun(input-n+2);\
			XMM4 = _mm_loadfun(input+n-2);\
			XMM7 = _mm_loadfun(input+n-4); /*out XMM7*/\
			XMM3 = _mm_shuffle_pd(XMM3, XMM4, _MM_SHUFFLE2(1, 0)); /*in XMM3*/\
			XMM4 = _mm_shuffle_pd(XMM4, XMM7, _MM_SHUFFLE2(1, 0));\
			XMM3 = _mm_add_pd(XMM3, XMM2);\
			XMM4 = _mm_add_pd(XMM4, XMM6);\
			XMM3 = _mm_mul_pd(XMM3, XMM1);\
			XMM4 = _mm_mul_pd(XMM4, XMM5);\
			XMM0 = _mm_add_pd(XMM0, XMM3);\
			XMM0 = _mm_add_pd(XMM0, XMM4)

#define SUM2(n, _mm_loadfun) \
			XMM1 = _mm_load_pd(coeff-n);\
			XMM5 = _mm_load_pd(coeff-n+2);\
			XMM2 = _mm_loadfun(input-n);\
			XMM6 = _mm_loadfun(input-n+2);\
			XMM4 = _mm_loadfun(input+n-2);\
			XMM3 = _mm_loadfun(input+n-4); /*out XMM3*/\
			XMM7 = _mm_shuffle_pd(XMM7, XMM4, _MM_SHUFFLE2(1, 0));/*in XMM7*/\
			XMM4 = _mm_shuffle_pd(XMM4, XMM3, _MM_SHUFFLE2(1, 0));\
			XMM7 = _mm_add_pd(XMM7, XMM2);\
			XMM4 = _mm_add_pd(XMM4, XMM6);\
			XMM7 = _mm_mul_pd(XMM7, XMM1);\
			XMM4 = _mm_mul_pd(XMM4, XMM5);\
			XMM0 = _mm_add_pd(XMM0, XMM7);\
			XMM0 = _mm_add_pd(XMM0, XMM4)



#define SUM01(n, _mm_loadfun) \
			XMM1 = _mm_load_pd(coeff-n);\
			XMM5 = _mm_load_pd(coeff-n+2);\
			\
			XMM2 = _mm_loadfun(input-n);\
			XMM6 = _mm_loadfun(input-n+2);\
			\
			XMM0 = _mm_load_sd(input+n);\
			XMM4 = _mm_loadfun(input+n-2);\
			XMM7 = _mm_loadfun(input+n-4); /*out XMM7*/\
			\
			XMM0 = _mm_shuffle_pd(XMM0, XMM4, _MM_SHUFFLE2(1, 0));\
			XMM4 = _mm_shuffle_pd(XMM4, XMM7, _MM_SHUFFLE2(1, 0));\
			XMM0 = _mm_add_pd(XMM0, XMM2);\
			XMM4 = _mm_add_pd(XMM4, XMM6);\
			XMM0 = _mm_mul_pd(XMM0, XMM1);\
			XMM4 = _mm_mul_pd(XMM4, XMM5);\
			XMM0 = _mm_add_pd(XMM0, XMM4)



#define ADD2BEGIN \
		XMM0 = _mm_load_pd(coef);\
		XMM1 = _mm_lddqu_pd(at);\
		XMM2 = _mm_load_pd(coef+2);\
		XMM3 = _mm_lddqu_pd(at+2);\
		\
		XMM1 = _mm_mul_pd(XMM1, XMM0);\
		XMM3 = _mm_mul_pd(XMM3, XMM2);\
		XMM3 = _mm_add_pd(XMM3, XMM1)

#define ADD2(n, XMM0, XMM1, XMM2, XMM3, XMM4) \
		XMM0 = _mm_load_pd(coef+n);\
		XMM1 = _mm_lddqu_pd(at+n);\
		XMM2 = _mm_load_pd(coef+n+2);\
		XMM3 = _mm_lddqu_pd(at+n+2);\
		\
		XMM1 = _mm_mul_pd(XMM1, XMM0);\
		XMM3 = _mm_mul_pd(XMM3, XMM2);\
		XMM3 = _mm_add_pd(XMM3, XMM4);\
		XMM3 = _mm_add_pd(XMM3, XMM1)

#define ADD3(n, XMM0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6) \
		XMM0 = _mm_load_pd(coef+n);\
		XMM1 = _mm_lddqu_pd(at+n);\
		XMM2 = _mm_load_pd(coef+n+2);\
		XMM3 = _mm_lddqu_pd(at+n+2);\
		XMM4 = _mm_load_pd(coef+n+4);\
		XMM5 = _mm_lddqu_pd(at+n+4);\
		\
		XMM1 = _mm_mul_pd(XMM1, XMM0);\
		XMM3 = _mm_mul_pd(XMM3, XMM2);\
		XMM5 = _mm_mul_pd(XMM5, XMM4);\
		\
		XMM3 = _mm_add_pd(XMM3, XMM1);\
		XMM5 = _mm_add_pd(XMM5, XMM6);\
		XMM5 = _mm_add_pd(XMM5, XMM3)

#define OUT(XMM) \
		XMM = _mm_hadd_pd(XMM, XMM);\
		_mm_store_sd(output+i, XMM)

#define ADD2_to3(n0) ADD2BEGIN
#define ADD2_3to4(n) ADD2(n, XMM0, XMM1, XMM2, XMM4, XMM3)
#define ADD2_4to3(n) ADD2(n, XMM0, XMM1, XMM2, XMM3, XMM4)
#define ADD3_3to6(n) ADD3(n, XMM0, XMM1, XMM2, XMM4, XMM5, XMM6, XMM3)
#define ADD3_4to6(n) ADD3(n, XMM0, XMM1, XMM2, XMM3, XMM5, XMM6, XMM4)
#define OUT_3 OUT(XMM3)
#define OUT_4 OUT(XMM4)
#define OUT_6 OUT(XMM6)


#define IADD1() \
	XMM2 = _mm_load_pd(coef);\
	XMM6 = _mm_load_pd(coef+2);\
	XMM0 = _mm_lddqu_pd(at);\
	XMM2 = _mm_mul_pd(XMM2, XMM1);\
	XMM6 = _mm_add_pd(XMM6, XMM2);\
	XMM0 = _mm_mul_pd(XMM0, XMM6)

#define IADD2(n) \
	XMM2 = _mm_load_pd(coef+2*n);\
	XMM3 = _mm_load_pd(coef+2*n+2);\
	XMM4 = _mm_load_pd(coef+2*n+4);\
	XMM5 = _mm_load_pd(coef+2*n+6);\
	XMM6 = _mm_lddqu_pd(at+n);\
	XMM7 = _mm_lddqu_pd(at+n+2);\
	\
	XMM2 = _mm_mul_pd(XMM2, XMM1);\
	XMM4 = _mm_mul_pd(XMM4, XMM1);\
	XMM3 = _mm_add_pd(XMM3, XMM2);\
	XMM5 = _mm_add_pd(XMM5, XMM4);\
	XMM6 = _mm_mul_pd(XMM6, XMM3);\
	XMM7 = _mm_mul_pd(XMM7, XMM5);\
	XMM0 = _mm_add_pd(XMM0, XMM6);\
	XMM0 = _mm_add_pd(XMM0, XMM7)


#define IADD3(n) \
	XMM3 = _mm_load_pd(coef+3*n);\
	XMM4 = _mm_load_pd(coef+3*n+2);\
	XMM5 = _mm_load_pd(coef+3*n+4);\
	XMM6 = _mm_lddqu_pd(at+n);\
	XMM3 = _mm_mul_pd(XMM3, XMM2);\
	XMM4 = _mm_mul_pd(XMM4, XMM1);\
	XMM5 = _mm_add_pd(XMM5, XMM3);\
	XMM5 = _mm_add_pd(XMM5, XMM4);\
	XMM6 = _mm_mul_pd(XMM6, XMM5);\
	XMM0 = _mm_add_pd(XMM0, XMM6)


