#pragma once
#include "Vector.hpp"
#include "Decomposition.hpp"
#include "SameBitCountConversion.hpp"
#include <type_traits>

#define LOG_PROCESS 0
#if (LOG_PROCESS == 1)
#include <iostream>
#endif

using namespace AppliedGeometry;

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
