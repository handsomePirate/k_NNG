#pragma once
#include "Vector.hpp"
#include <type_traits>

// Templates to convert types of the same bit count.
// {
template<typename T>
struct SameBitCount {};

#define SAME_BIT_CONSTRUCT(T1, T2) \
template<> struct SameBitCount<T1> { using type = T2; }; \
template<> struct SameBitCount<T2> { using type = T1; }; \

SAME_BIT_CONSTRUCT(float, int32_t)
SAME_BIT_CONSTRUCT(double, int64_t)

#define eqtype(T) typename SameBitCount<T>::type

template<typename T>
union SameBitConverter
{
	T x1;
	eqtype(T) x2;
};
// }

using namespace AppliedGeometry;

// This class is capable of comparing two n-dimensional points & determining their z-order.
class Morton
{
public:
	// This is the main comparison function, it finds the highest most significant bit of all
	// dimensions and determines the order of two points based on the comparison of the elements
	// in that dimension.
	template<int N, typename T>
	static bool Compare(const VectorN<N, T>& v1, const VectorN<N, T>& v2)
	{
		int32_t x = 0;
		int32_t d = 0;
		int32_t y;
		for (int32_t n = 0; n < N; ++n)
		{
			y = XORMSB<T>(v1.GetComponent(n), v2.GetComponent(n));
			if (x < y)
			{
				x = y;
				d = n;
			}
		}
		return v1.GetComponent(d) < v2.GetComponent(d);
	}
	// This function allows to plug the Morton class into a sort.
	template<int N, typename T>
	bool operator()(const VectorN<N, T>& v1, const VectorN<N, T>& v2) const
	{
		return Compare<N, T>(v1, v2);
	}

private:
	// For floating point numbers, we need to take a different approach.
	template<typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
	static int32_t XORMSB(const T& a, const T& b)
	{
		int32_t exponentA;
		eqtype(T) mantissaA = GetFloatComponents(a, exponentA);

		int32_t exponentB;
		eqtype(T) mantissaB = GetFloatComponents(b, exponentB);

		if (exponentA == exponentB)
		{
			if (mantissaA == mantissaB)
			{
				return 0;
			}
			int32_t z = MSB<eqtype(T)>(mantissaA, mantissaB);
			return exponentA - 24 + z;
		}
		else
		{
			return (exponentB < exponentA) ? exponentA : exponentB;
		}
	}
	
	// For integers, we can compute what we need directly.
	template<typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
	static int32_t XORMSB(const T& a, const T& b)
	{
		return MSB<T>(a, b);
	}

	// Computes the most significant bit position.
	template<typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
	static int32_t MSB(const T& a, const T& b)
	{
		T x = a ^ b;
		int32_t msb = 0;
		while (x != 0) 
		{
			x >>= 1;
			msb++;
		}
		return msb;
	}

	// This function offers a view of the bits in 'x' through the type 'T' is paired with.
	template<typename T>
	static eqtype(T) ConvertToSameBitCount(const T& x)
	{
		SameBitConverter<T> c;
		c.x1 = x;
		return c.x2;
	}

	template<typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
	static eqtype(T) GetFloatComponents(const T& x, int32_t& e)
	{
		eqtype(T) i = ConvertToSameBitCount(x);
		e = ((0x7F800000 & i) >> 23);
		e = e == 0 ? 0 : e - 126;
		return 0x007FFFFF & i;
	}
	// There are different versions of frexp functions for floats and doubles,
	// this requires an extra level of abstraction.
	static float Components(float x, int32_t& e)
	{
		return frexpf(x, &e);
	}
	static double Components(double x, int32_t& e)
	{
		return frexp(x, &e);
	}
}; // 1 8 23
// mant: 0x007FFFFF
// expo: 0x7F800000
// sign: 0x80000000
// s&e:  0xFF800000
