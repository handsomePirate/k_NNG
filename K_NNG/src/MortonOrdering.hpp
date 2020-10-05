#pragma once
#include "Vector.hpp"
#include <type_traits>

// 0000 0000 0000 0000 0000 0000 0001 0111
// 0x007FFFFF
typedef union {
	float x;
	struct {
		uint64_t mantissa : 23;
		uint32_t exponent : 8;
		uint32_t sign : 1;
	} parts;
} float_cast;

typedef union {
	double x;
	struct {
		uint64_t mantissa : 52;
		uint32_t exponent : 11;
		uint32_t sign : 1;
	} parts;
} double_cast;

using namespace AppliedGeometry;

class Morton
{
public:
	template<int N, typename T>
	static bool Compare(const VectorN<N, T>& v1, const VectorN<N, T>& v2)
	{
		uint32_t x = 0;
		uint32_t d = 0;
		uint32_t y;
		for (uint32_t n = 0; n < N; ++n)
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
	template<int N, typename T>
	bool operator()(const VectorN<N, T>& v1, const VectorN<N, T>& v2)
	{
		return Compare<N, T>(v1, v2);
	}

private:
	template<typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
	static uint32_t XORMSB(const T& a, const T& b)
	{
		uint32_t exponentA;
		uint64_t mantissaA;
		Components<T>(a, mantissaA, exponentA);

		uint32_t exponentB;
		uint64_t mantissaB;
		Components<T>(b, mantissaB, exponentB);

		if (exponentA == exponentB)
		{
			uint32_t z = MSB<uint64_t>(mantissaA, mantissaB);
			return exponentA - z;
		}
		else
		{
			return (exponentB < exponentA) ? exponentA : exponentB;
		}
	}
	
	template<typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
	static uint32_t XORMSB(const T& a, const T& b)
	{
		return MSB<T>(a, b);
	}

	
	template<typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
	static uint32_t MSB(const T& a, const T& b)
	{
		T x = a ^ b;
		uint32_t msb = 0;
		while (x != 0) 
		{
			x >>= 1;
			msb++;
		}
		return msb;
	}

	template<typename T, std::enable_if_t<std::is_same<T, float>::value, int> = 0>
	static void Components(T x, uint64_t& m, uint32_t& e)
	{
		m = (uint32_t&)frexpf(x, &e);
	}
	template<typename T, std::enable_if_t<std::is_same<T, double>::value, int> = 0>
	static void Components(T x, uint64_t& m, uint32_t& e)
	{
		m = (uint64_t&)frexp(x, &e);
	}
};