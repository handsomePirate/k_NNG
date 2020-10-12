#pragma once
#include <stdint.h>

// This allows for a floating point number's decomposition - provides masks, etc.

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
