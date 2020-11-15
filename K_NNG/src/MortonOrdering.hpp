#pragma once
#include "Vector.hpp"
#include "Decomposition.hpp"
#include "SameBitCountConversion.hpp"
#include <type_traits>
#include <iostream>

#define LOG_PROCESS 0
#if (LOG_PROCESS == 1)
#include <iostream>
#endif

using namespace AppliedGeometry;

inline int lzcount32(uint32_t x)
{
	int r = 0;

	while (x >>= 1)
		r++;
	return r;
}

inline int lzcount64(uint64_t x)
{
	int r = 0;

	while (x >>= 1)
		r++;
	return r;
}

typedef union { float f; int32_t i; uint32_t u; } union32_t;
typedef union { double f; int64_t i; uint64_t u; } union64_t;

inline uint32_t ieee_exponent(float f)
{
	union32_t u;
	u.f = f;
	return (u.i & 0x7f800000) >> 23;
}

inline uint32_t ieee_mantissa(float f)
{
	union32_t u;
	u.f = f;
	return (u.i & 0x007FFFFF);
}

inline uint64_t ieee_exponent(double f)
{
	union64_t u;
	u.f = f;
	return (u.i & 0x7ff0000000000000ll) >> 52;
}

inline uint64_t ieee_mantissa(double f)
{
	union64_t u;
	u.f = f;
	return (u.i & 0x000fffffffffffffll);
}

void WriteBinFloat(float x, std::ostream& os)
{
	unsigned int* px = (unsigned int*)&x;
	os << '(';
	for (unsigned int b = 0x80000000; b > 0; b >>= 1)
	{
		if (*px & b)
		{
			os << '1';
		}
		else
		{
			os << '0';
		}
		if (b == 0x80000000 || b == 0x00800000)
		{
			os << ' ';
		}
	}
	os << " -> ";
	int exp;
	frexpf(x, &exp);
	os << "exp=" << exp;
	os << ')';
}

template<int ND, class FloatType>
class FloatMortonLess
{
public:
	bool operator() (const VectorN<ND, FloatType>& v1, const VectorN<ND, FloatType>& v2) const
	{
		int x(0), dim(0);

		for (int j = 0; j < ND; ++j) 
		{
			int y = xormsb(v1.GetComponent(j), v2.GetComponent(j));

			if (x < y) 
			{
				x = y;
				dim = j;
			}
		}

		return v1.GetComponent(dim) < v2.GetComponent(dim);
	}

private:
	int xormsb(FloatType a, FloatType b) const
	{
		int x = ieee_exponent(a);
		int y = ieee_exponent(b);

		if (x == y)
		{
			int z = msdb(ieee_mantissa(a), ieee_mantissa(b));
			return x - z;
		}

		return (y < x) ? x : y;
	}

	int msdb(uint32_t a, uint32_t b) const 
	{
		return 32 - lzcount32(a ^ b);
	}

	int msdb(uint64_t a, uint64_t b) const 
	{
		return 64 - lzcount64(a ^ b);
	}
};

// The idea here is to enable comparison of two n-dimensional vectors by comparing their coordinates
// in the dimension that has the highest most significant bit after XORing those two values.
// For integers, this is straight forward, real numbers have to be considered in a fixed point representation.
// The format of a floating point number is (-1)^s * m * 2^e, where
// s = sign,
// m = mantissa,
// e = exponent.
// Thus, the exponent directly controls the position of the most significant bit in fixed point representation.
// If the exponents are not equal, the greater one is the resulting MSB of the XOR. If they are the same, 
// it is necessary to XOR the mantissa and add the resulting bit shift to one of the exponents.
// Ordered in this way, the vectors that are nearby in a linear memory arrangement are also close in the
// n-dimensional space.
// Therefore, when working with k nearest neighbours, we want precisely this locality, so that different threads access
// memory in fewer places, which makes the caches very effective.

// This class is capable of comparing two n-dimensional points & determining their z-order.
template<int N, typename T>
class Morton
{
public:
	// A set of minimum values for each dimension so that the ordering algorithm works only
	// with positive numbers.
	Morton(T* mins, size_t minsSize)
	{
		assert(N == minsSize);
		for (int n = 0; n < N; ++n)
		{
			mins_[n] = mins[n];
		}
	}

	Morton()
	{
		for (int n = 0; n < N; ++n)
		{
			mins_[n] = 0;
		}
	}

	// This function allows to plug the Morton class into a sort.
	bool operator()(const VectorN<N, T>& v1, const VectorN<N, T>& v2) const
	{
		return Compare(v1, v2);
	}

private:
	// The minimum values of each dimension (needs to be passed into this structure).
	// They are used to avoid working with negative numbers.
	std::array<T, N> mins_;

	// This is the main comparison function, it finds the highest most significant bit of all
	// dimensions and determines the order of two points based on the comparison of the elements
	// in that dimension.
	bool Compare(const VectorN<N, T>& v1, const VectorN<N, T>& v2) const
	{
		int32_t x = INT32_MIN;
		int32_t d = 0;
		int32_t y;
		for (int32_t n = 0; n < N; ++n)
		{
			//float x1 = v1.GetComponent(n) - mins_[n];
			//float x2 = v2.GetComponent(n) - mins_[n];
			//WriteBinFloat(x1, std::cout);
			//std::cout << ',';
			//WriteBinFloat(x2, std::cout);
			//std::cout << std::endl;
			y = XORMSB<T>(v1.GetComponent(n) - mins_[n], v2.GetComponent(n) - mins_[n]);
			if (x < y)
			{
				x = y;
				d = n;
			}
	}
#if (LOG_PROCESS == 1)
		std::cout << v1.ToString() << " x " << v2.ToString() << " => " << (v1.GetComponent(d) < v2.GetComponent(d)) << std::endl;
#endif
		//std::cout << std::endl;
		return v1.GetComponent(d) < v2.GetComponent(d);
}

	// For floating point numbers, we need to take a different approach.
	template<typename Ty, std::enable_if_t<std::is_floating_point<Ty>::value, int> = 0>
	static int32_t XORMSB(const Ty& a, const Ty& b)
	{
		if (a == b)
		{
			return INT32_MIN;
		}

		int32_t exponentA;
		eqtype(Ty) mantissaA = GetFloatComponents(a, exponentA);

		int32_t exponentB;
		eqtype(Ty) mantissaB = GetFloatComponents(b, exponentB);

		const int32_t exponentComplement = Decomposition<Ty>::ExponentComplement();
		if (exponentA == exponentB)
		{
			int32_t z = MSB<eqtype(Ty)>(mantissaA ^ mantissaB);
			return exponentA - exponentComplement - Decomposition<Ty>::MantissaBitCount() + z;
		}
		else
		{
			if (a == 0)
			{
				return exponentB - exponentComplement + 1;
			}
			if (b == 0)
			{
				return exponentA - exponentComplement + 1;
			}
			return (exponentB < exponentA) ? exponentA - exponentComplement + 1 : exponentB - exponentComplement + 1;
		}
	}
	
	// For integers, we can compute what we need directly.
	template<typename Ty, std::enable_if_t<std::is_integral<Ty>::value, int> = 0>
	static int32_t XORMSB(const Ty& a, const Ty& b)
	{
		return MSB<Ty>(a ^ b);
	}

	// Computes the most significant bit position.
	template<typename Ty, std::enable_if_t<std::is_integral<Ty>::value, int> = 0>
	static int32_t MSB(Ty&& x)
	{
		int32_t msb = 0;
		while (x != 0) 
		{
			x >>= 1;
			msb++;
		}
		return msb;
	}

	// This function offers a view of the bits in 'x' through the type 'Ty' is paired with.
	template<typename Ty>
	static eqtype(Ty) ConvertToSameBitCount(const Ty& x)
	{
		SameBitConvertor<Ty> c;
		c.x1 = x;
		return c.x2;
	}

	template<typename Ty, std::enable_if_t<std::is_floating_point<Ty>::value, int> = 0>
	static eqtype(Ty) GetFloatComponents(const Ty& x, int32_t& e)
	{
		eqtype(Ty) i = ConvertToSameBitCount(x);
		e = ((Decomposition<Ty>::ExponentMask() & i) >> Decomposition<Ty>::MantissaBitCount());
		return Decomposition<Ty>::MantissaMask() & i; // TODO: float -> double
	}
};
