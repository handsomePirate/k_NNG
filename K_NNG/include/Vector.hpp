#pragma once
#include <type_traits>
#include <string>
#include <sstream>
#include <cassert>
#include <array>

namespace AppliedGeometry
{
	class PointDistance;

	template <typename NumericType> struct Vector2;
	template <typename NumericType> struct Vector3;
	template <typename NumericType> struct Vector4;

	/// A general N-dimensional point.
	template <size_t N, typename NumericType,
		std::enable_if_t<std::is_arithmetic<NumericType>::value, int> = 0>
	struct VectorN
	{
		/// Initializer list constructor.
		template <typename ... List>
		VectorN(List&& ... components)
			: components_{ { std::forward<List>(components)... } } {}
		/// Default constructor.
		VectorN() : components_{} {}

		/// Copy & Move constructors & operators.
		VectorN(VectorN&& other) noexcept { std::swap(components_, other.components_); }
		VectorN(const VectorN& other) { components_ = other.components_; }
		VectorN& operator=(VectorN&& other) 
		{ 
			std::swap(components_, other.components_);
			return *this;
		}
		VectorN& operator=(const VectorN& other)
		{
			components_ = other.components_;
			return *this;
		}

		/// Backwards conversions.
		VectorN(Vector2<NumericType>&& other) 
			: components_(std::move(other.components_)) 
		{
			static_assert(N == 2, "Cannot convert a 2-dimensional vector to other dimensions implicitly.");
		}
		VectorN(const Vector2<NumericType>& other)
			: components_(other.components_)
		{
			static_assert(N == 2, "Cannot convert a 2-dimensional vector to other dimensions implicitly.");
		}

		VectorN(Vector3<NumericType>&& other)
			: components_(std::move(other.components_))
		{
			static_assert(N == 3, "Cannot convert a 3-dimensional vector to other dimensions implicitly.");
		}
		VectorN(const Vector3<NumericType>& other)
			: components_(other.components_)
		{
			static_assert(N == 3, "Cannot convert a 3-dimensional vector to other dimensions implicitly.");
		}

		VectorN(Vector4<NumericType>&& other)
			: components_(std::move(other.components_))
		{
			static_assert(N == 4, "Cannot convert a 4-dimensional vector to other dimensions implicitly.");
		}
		VectorN(const Vector4<NumericType>& other)
			: components_(other.components_)
		{
			static_assert(N == 4, "Cannot convert a 4-dimensional vector to other dimensions implicitly.");
		}
		
		/// Modifiable component access.
		NumericType& operator[](const size_t i)
		{
			assert(i < N);
			return components_[i];
		}

		/// Const component access.
		const NumericType& GetComponent(const size_t i) const
		{
			assert(i < N);
			return components_[i];
		}

		/// Distance to another point.
		NumericType Distance(const VectorN<N, NumericType>& other) const
		{
			return PointDistance::Get(*this, other);
		}

		/// Distance to another point squared.
		NumericType DistanceSqr(const VectorN<N, NumericType>& other) const
		{
			return PointDistance::GetSqr(*this, other);
		}

		/// Converts the point to a different format.
		template <typename NewNumericType,
			std::enable_if_t<!std::is_same_v<NewNumericType, NumericType>, int> = 0>
		VectorN<N, NewNumericType> To() const
		{
			VectorN<N, NewNumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = NewNumericType(components_[i]);
			}
			return result;
		}

		/// Implements a default * operator for vectors.
		VectorN<N, NumericType> operator*(const VectorN<N, NumericType>& other) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] * other.GetComponent(i);
			}
			return result;
		}

		/// Implements a default + operator for vectors.
		VectorN<N, NumericType> operator+(const VectorN<N, NumericType>& other) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] + other.GetComponent(i);
			}
			return result;
		}

		/// Implements a default - operator for vectors.
		VectorN<N, NumericType> operator-(const VectorN<N, NumericType>& other) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] - other.GetComponent(i);
			}
			return result;
		}

		/// Implements a default / operator for vectors.
		VectorN<N, NumericType> operator/(const VectorN<N, NumericType>& other) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] / other.GetComponent(i);
			}
			return result;
		}

		/// Implements a default * operator for scalars.
		VectorN<N, NumericType> operator*(const NumericType& scalar) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] * scalar;
			}
			return result;
		}

		/// Implements a default + operator for scalars.
		VectorN<N, NumericType> operator+(const NumericType& scalar) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] + scalar;
			}
			return result;
		}

		/// Implements a default - operator for scalars.
		VectorN<N, NumericType> operator-(const NumericType& scalar) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result[i] = components_[i] - scalar;
			}
			return result;
		}

		/// Implements a default / operator for scalars.
		VectorN<N, NumericType> operator/(const NumericType& scalar) const
		{
			VectorN<N, NumericType> result{};
			for (size_t i = 0; i < N; ++i)
			{
				result.Component(i) = components_[i] / scalar;
			}
			return result;
		}

		/// Implements a default *= operator for vectors.
		VectorN<N, NumericType>& operator*=(const VectorN<N, NumericType>& other)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] *= other.GetComponent(i);
			}
			return *this;
		}

		/// Implements a default += operator for vectors.
		VectorN<N, NumericType>& operator+=(const VectorN<N, NumericType>& other)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] += other.GetComponent(i);
			}
			return *this;
		}

		/// Implements a default -= operator for vectors.
		VectorN<N, NumericType>& operator-=(const VectorN<N, NumericType>& other)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] -= other.GetComponent(i);
			}
			return *this;
		}

		/// Implements a default /= operator for vectors.
		VectorN<N, NumericType>& operator/=(const VectorN<N, NumericType>& other)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] /= other.GetComponent(i);
			}
			return *this;
		}

		/// Implements a default *= operator for scalars.
		VectorN<N, NumericType>& operator*=(const NumericType& scalar)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] *= scalar;
			}
			return *this;
		}

		/// Implements a default += operator for scalars.
		VectorN<N, NumericType>& operator+=(const NumericType& scalar)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] += scalar;
			}
			return *this;
		}

		/// Implements a default -= operator for scalars.
		VectorN<N, NumericType>& operator-=(const NumericType& scalar)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] -= scalar;
			}
			return *this;
		}

		/// Implements a default /= operator for scalars.
		VectorN<N, NumericType>& operator/=(const NumericType& scalar)
		{
			for (size_t i = 0; i < N; ++i)
			{
				components_[i] /= scalar;
			}
			return *this;
		}

		/// Computes the size of the vector.
		NumericType Size() const
		{
			return PointDistance::Get(*this, VectorN{});
		}

		/// Normalizes the vector.
		void Normalize()
		{
			*this /= Size();
		}

		/// Normalizes the vector.
		VectorN Normalized()
		{
			return *this / Size();
		}

		/// Creates a hash from all the components of the vector.
		size_t Hash() const
		{
			size_t hash = 0;
			// Allocate a memory, that the bit shifts will happen in.
			size_t componentMemory = 0;
			// Have a view of the memory in the Vector data type.
			NumericType* componentView = (NumericType*)&componentMemory;
			// Hash each component.
			for (size_t i = 0; i < N; ++i)
			{
				const int bitMove = (int)sizeof(NumericType) + 1;
				*componentView = components_[i];
				// Compute and use a rotating shift.
				const int shift = (bitMove * i) % sizeof(unsigned long long);
				auto moved = componentMemory << (shift);
				const size_t bitCount = sizeof(unsigned long long) * 8;
				auto movedComplement = 
					componentMemory >> (bitCount - shift);
				hash ^= moved | movedComplement;
				componentMemory = 0;
			}
			return hash;
		}

		/// An automatic conversion to string.
		std::string ToString() const
		{
			std::stringstream ss;
			ss << '(';
			size_t element = 0;
			do
			{
				ss << components_[element];
			} while (++element < N && (ss << ", "));
			ss << ')';
			return ss.str();
		}
	protected:
		/// The components of the N-dimensional point.
		std::array<NumericType, N> components_;
	};

	/// Implements a default * operator for scalars.
	template <size_t N, typename NumericType,
		std::enable_if_t<std::is_arithmetic<NumericType>::value, int> = 0>
	VectorN<N, NumericType> operator*(const NumericType& scalar, const VectorN<N, NumericType>& vector)
	{
		VectorN<N, NumericType> result{};
		for (size_t i = 0; i < N; ++i)
		{
			result[i] = scalar * vector.GetComponent(i);
		}
		return result;
	}

	/// Implements a default + operator for scalars.
	template <size_t N, typename NumericType,
		std::enable_if_t<std::is_arithmetic<NumericType>::value, int> = 0>
	VectorN<N, NumericType> operator+(const NumericType& scalar, const VectorN<N, NumericType>& vector)
	{
		VectorN<N, NumericType> result{};
		for (size_t i = 0; i < N; ++i)
		{
			result[i] = scalar + vector.GetComponent(i);
		}
		return result;
	}

	/// A 2-dimensional point.
	template <typename NumericType>
	struct Vector2 : public VectorN<2, NumericType>
	{
		/// Initialization constructor.
		Vector2(const NumericType& x, const NumericType& y) : VectorN<2, NumericType>(x, y) {}
		/// Default constructor.
		Vector2() : VectorN<2, NumericType>{} {}

		/// General vector constructors.
		Vector2(VectorN<2, NumericType>&& other) : VectorN<2, NumericType>(std::move(other)) {}
		Vector2(const VectorN<2, NumericType>& other) : VectorN<2, NumericType>(other) {}

		/// Modifiable component access.
		NumericType& X() { return VectorN<2, NumericType>::operator[](0); }
		NumericType& Y() { return VectorN<2, NumericType>::operator[](1); }

		/// Const component access.
		const NumericType& GetX() const { return VectorN<2, NumericType>::GetComponent(0); }
		const NumericType& GetY() const { return VectorN<2, NumericType>::GetComponent(1); }

		/// Converts the point to a different format.
		template <typename NewNumericType,
			std::enable_if_t<!std::is_same_v<NewNumericType, NumericType>, int> = 0>
		Vector2<NewNumericType> To() const
		{
			return Vector2<NewNumericType>(NewNumericType(GetX()), NewNumericType(GetY()));
		}

		/// Implements a * operator for 2-dimensional vectors.
		Vector2<NumericType> operator*(const Vector2<NumericType>& other) const
		{
			return Vector2<NumericType>(GetX() * other.GetX(), GetY() * other.GetY());
		}

		NumericType Cross(const Vector2<NumericType>& other) const
		{
			return GetX() * other.GetY() - GetY() * other.GetX();
		}
	};

	/// A 3-dimensional point.
	template <typename NumericType>
	struct Vector3 : public VectorN<3, NumericType>
	{
		/// Initialization constructor.
		Vector3(const NumericType& x, const NumericType& y, const NumericType& z) 
			: VectorN<3, NumericType>(x, y, z) {}
		/// Default constructor.
		Vector3() : VectorN<3, NumericType>{} {}

		/// General vector constructors.
		Vector3(VectorN<3, NumericType>&& other) : VectorN<3, NumericType>(std::move(other)) {}
		Vector3(const VectorN<3, NumericType>& other) : VectorN<3, NumericType>(other) {}

		/// Modifiable component access.
		NumericType& X() { return VectorN<3, NumericType>::operator[](0); }
		NumericType& Y() { return VectorN<3, NumericType>::operator[](1); }
		NumericType& Z() { return VectorN<3, NumericType>::operator[](2); }

		/// Const component access.
		const NumericType& GetX() const { return VectorN<3, NumericType>::GetComponent(0); }
		const NumericType& GetY() const { return VectorN<3, NumericType>::GetComponent(1); }
		const NumericType& GetZ() const { return VectorN<3, NumericType>::GetComponent(2); }

		/// Creates a 2-dimensional Point from the x and y components.
		template <typename NewNumericType = NumericType>
		Vector2<NewNumericType> XY() const 
		{ 
			return Vector2<NewNumericType>(NewNumericType(GetX()), NewNumericType(GetY())); 
		}

		/// Creates a 2-dimensional Point from the x and z components.
		template <typename NewNumericType = NumericType>
		Vector2<NewNumericType> XZ() const 
		{ 
			return Vector2<NewNumericType>(NewNumericType(GetX()), NewNumericType(GetZ()));
		}

		/// Creates a 2-dimensional Point from the y and z components.
		template <typename NewNumericType = NumericType>
		Vector2<NewNumericType> YZ() const 
		{ 
			return Vector2<NewNumericType>(NewNumericType(GetY()), NewNumericType(GetZ()));
		}

		/// Converts the point to a different format.
		template <typename NewNumericType,
			std::enable_if_t<!std::is_same_v<NewNumericType, NumericType>, int> = 0>
		Vector3<NewNumericType> To() const
		{
			return Vector3<NewNumericType>(
				NewNumericType(GetX()), NewNumericType(GetY()), NewNumericType(GetZ()));
		}

		/// Implements a * operator for 3-dimensional vectors.
		Vector3<NumericType> operator*(const Vector3<NumericType>& other) const
		{
			return Vector3<NumericType>(GetX() * other.GetX(), GetY() * other.GetY(), GetZ() * other.GetZ());
		}
	};

	/// A 4-dimensional point.
	template <typename NumericType>
	struct Vector4 : public VectorN<4, NumericType>
	{
		/// Initialization constructor.
		Vector4(const NumericType& x, const NumericType& y, const NumericType& z, const NumericType& w)
			: VectorN<4, NumericType>({ x, y, z }) {}
		/// Default constructor.
		Vector4() : VectorN<4, NumericType>{} {}

		/// General vector constructors.
		Vector4(VectorN<4, NumericType>&& other) : VectorN<4, NumericType>(std::move(other)) {}
		Vector4(const VectorN<4, NumericType>& other) : VectorN<4, NumericType>(other) {}

		/// Modifiable component access.
		NumericType& X() { return VectorN<4, NumericType>::operator[](0); }
		NumericType& Y() { return VectorN<4, NumericType>::operator[](1); }
		NumericType& Z() { return VectorN<4, NumericType>::operator[](2); }
		NumericType& W() { return VectorN<4, NumericType>::operator[](3); }

		/// Const component access.
		const NumericType& GetX() const { return VectorN<4, NumericType>::GetComponent(0); }
		const NumericType& GetY() const { return VectorN<4, NumericType>::GetComponent(1); }
		const NumericType& GetZ() const { return VectorN<4, NumericType>::GetComponent(2); }
		const NumericType& GetW() const { return VectorN<4, NumericType>::GetComponent(3); }

		/// Creates a 2-dimensional Point from the x and y components.
		template <typename NewNumericType = NumericType>
		Vector2<NewNumericType> XY() const
		{
			return Vector2<NewNumericType>(NewNumericType(GetX()), NewNumericType(GetY()));
		}

		/// Creates a 2-dimensional Point from the x and z components.
		template <typename NewNumericType = NumericType>
		Vector2<NewNumericType> XZ() const
		{
			return Vector2<NewNumericType>(NewNumericType(GetX()), NewNumericType(GetZ()));
		}

		/// Creates a 2-dimensional Point from the y and z components.
		template <typename NewNumericType = NumericType>
		Vector2<NewNumericType> YZ() const
		{
			return Vector2<NewNumericType>(NewNumericType(GetY()), NewNumericType(GetZ()));
		}

		/// Creates a 3-dimensional Point from the x, y and z components.
		template <typename NewNumericType = NumericType>
		Vector3<NewNumericType> XYZ() const
		{
			return Vector3<NewNumericType>(
				NewNumericType(GetX()), NewNumericType(GetY()), NewNumericType(GetZ()));
		}

		/// Converts the point to a different format.
		template <typename NewNumericType,
			std::enable_if_t<!std::is_same_v<NewNumericType, NumericType>, int> = 0>
		Vector4<NewNumericType> To() const
		{
			return Vector4<NewNumericType>(
				NewNumericType(GetX()), NewNumericType(GetY()), NewNumericType(GetZ()), 
				NewNumericType(GetW()));
		}

		/// Implements a * operator for 4-dimensional vectors.
		Vector4<NumericType> operator*(const Vector4<NumericType>& other) const
		{
			return Vector4<NumericType>(GetX() * other.GetX(), GetY() * other.GetY(), 
				GetZ() * other.GetZ(), GetW() * other.GetW());
		}
	};

	/// Define 2-dimensional point types.
	typedef Vector2<int> Vector2I;
	typedef Vector2<float> Vector2F;
	typedef Vector2<double> Vector2D;

	/// Define 3-dimensional point types.
	typedef Vector3<int> Vector3I;
	typedef Vector3<float> Vector3F;
	typedef Vector3<double> Vector3D;

	/// Define 4-dimensional point types.
	typedef Vector4<int> Vector4I;
	typedef Vector4<float> Vector4F;
	typedef Vector4<double> Vector4D;

	/// This class statically computes either the distance or the distance squared.
	class PointDistance
	{
	public:
		/// Compute the distance between the two points squared.
		template <size_t N, typename NumericType>
		static NumericType GetSqr(
			const VectorN<N, NumericType>& point1, const VectorN<N, NumericType>& point2)
		{
			NumericType sum{};
			for (size_t i = 0; i < N; ++i)
			{
				NumericType sub = point1.GetComponent(i) - point2.GetComponent(i);
				sum += sub * sub;
			}
			return sum;
		}

		/// Compute the distance between the two points.
		template <size_t N, typename NumericType>
		static NumericType Get(
			const VectorN<N, NumericType>& point1, const VectorN<N, NumericType>& point2)
		{
			return NumericType(sqrt(GetSqr(point1, point2)));
		}
	};

	template <size_t N, typename NumericType>
	bool operator==(const VectorN<N, NumericType>& left, const VectorN<N, NumericType>& right)
	{
		for (size_t i = 0; i < N; ++i)
		{
			if (left.GetComponent(i) != right.GetComponent(i))
				return false;
		}
		return true;
	}
}

namespace std
{
	using namespace AppliedGeometry;

	template <size_t N, typename NumericType>
	struct hash<VectorN<N, NumericType>>
	{
		size_t operator()(const VectorN<N, NumericType>& vec) const noexcept
		{
			return vec.Hash();
		}
	};

	template <typename NumericType>
	struct hash<Vector2<NumericType>>
	{
		size_t operator()(const Vector2<NumericType>& vec) const noexcept
		{
			return vec.Hash();
		}
	};

	template <typename NumericType>
	struct hash<Vector3<NumericType>>
	{
		size_t operator()(const Vector3<NumericType>& vec) const noexcept
		{
			return vec.Hash();
		}
	};

	template <typename NumericType>
	struct hash<Vector4<NumericType>>
	{
		size_t operator()(const Vector4<NumericType>& vec) const noexcept
		{
			return vec.Hash();
		}
	};
}
