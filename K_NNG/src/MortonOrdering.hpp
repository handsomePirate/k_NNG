#pragma once
#include "Vector.hpp"
#include <type_traits>


#define LOG_PROCESS 0
#if (LOG_PROCESS == 1)
#include <iostream>
#endif

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

template<typename T>
struct Decomposition {};
template<>
struct Decomposition<float> 
{
	// mant: 0x007FFFFF
	// expo: 0x7F800000
	// sign: 0x80000000
	// s&e:  0xFF800000

	static inline int32_t MantissaMask()
	{
		return 0x007FFFFF;
	}
	static inline int32_t MantissaBitCount()
	{
		return 23;
	}
	static inline int32_t ExponentMask()
	{
		return 0x7F800000;
	}
	static inline int32_t ExponentComplement()
	{
		return 127;
	}
};
template<>
struct Decomposition<double> 
{
	// mant: 0x000FFFFFFFFFFFFF
	// expo: 0x7FF0000000000000
	// sign: 0x8000000000000000
	// s&e:  0xFFF0000000000000

	static inline int64_t MantissaMask()
	{
		return 0x000FFFFFFFFFFFFF;
	}
	static inline int32_t MantissaBitCount()
	{
		return 52;
	}
	static inline int64_t ExponentMask()
	{
		return 0x7FF0000000000000;
	}
	static inline int32_t ExponentComplement()
	{
		return 1024;
	}
};

using namespace AppliedGeometry;

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
		int32_t exponentA;
		eqtype(Ty) mantissaA = GetFloatComponents(a, exponentA);

		int32_t exponentB;
		eqtype(Ty) mantissaB = GetFloatComponents(b, exponentB);

		const int32_t exponentComplement = Decomposition<Ty>::ExponentComplement();
		if (exponentA == exponentB)
		{
			if (mantissaA == mantissaB)
			{
				return INT32_MIN;
			}
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

	// This function offers a view of the bits in 'x' through the type 'T' is paired with.
	template<typename Ty>
	static eqtype(Ty) ConvertToSameBitCount(const Ty& x)
	{
		SameBitConverter<Ty> c;
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
};
