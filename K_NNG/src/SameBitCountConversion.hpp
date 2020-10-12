#pragma once
#include <stdint.h>

// Templates to convert types of the same bit count.

template<typename T>
struct SameBitCount {};

#define SAME_BIT_CONSTRUCT(T1, T2) \
template<> struct SameBitCount<T1> { using type = T2; }; \
template<> struct SameBitCount<T2> { using type = T1; }; \

SAME_BIT_CONSTRUCT(float, int32_t)
SAME_BIT_CONSTRUCT(double, int64_t)

#define eqtype(T) typename SameBitCount<T>::type

template<typename T>
union SameBitConvertor
{
	T x1;
	eqtype(T) x2;
};
