#ifndef _dc7a30fd_6048_4dfc_a551_efae434269f0
#define _dc7a30fd_6048_4dfc_a551_efae434269f0

#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

namespace sycomore
{

/// @addtogroup Arrays
/// @{

/// @brief Real numbers
using Real = double;

/// @brief Complex numbers
using Complex = std::complex<Real>;

/// @brief 2D Vector of arithmetic types
template <typename T>
using Vector2 = xt::xtensor_fixed<T, xt::xshape<2>>;

/// @brief 3D Vector of arithmetic types
template <typename T>
using Vector3 = xt::xtensor_fixed<T, xt::xshape<3>>;

/// @brief 4D Vector of arithmetic types
template <typename T>
using Vector4 = xt::xtensor_fixed<T, xt::xshape<4>>;

/// @brief 2x2 matrix of arithmetic types
template <typename T>
using Matrix2x2 = xt::xtensor_fixed<T, xt::xshape<2, 2>>;

/// @brief 3x3 matrix of arithmetic types
template <typename T>
using Matrix3x3 = xt::xtensor_fixed<T, xt::xshape<3, 3>>;

/// @brief 4x4 matrix of arithmetic types
template <typename T>
using Matrix4x4 = xt::xtensor_fixed<T, xt::xshape<4, 4>>;

/// @brief Static-dimension array of real numbers
template<std::size_t N, xt::layout_type L=XTENSOR_DEFAULT_LAYOUT>
using TensorR = xt::xtensor<Real, N, L>;

/// @brief Static-dimension array of complex numbers
template<std::size_t N, xt::layout_type L=XTENSOR_DEFAULT_LAYOUT>
using TensorC = xt::xtensor<Complex, N, L>;

/// @brief Dynamic-dimension array of real numbers
using ArrayR = xt::xarray<Real>;

/// @brief Dynamic-dimension array of complex numbers
using ArrayC = xt::xarray<Complex>;

/// @brief 2D vector of real numbers
using Vector2R = Vector2<Real>;

/// @brief 3D vector of real numbers
using Vector3R = Vector3<Real>;

/// @brief 4D vector of real numbers
using Vector4R = Vector4<Real>;

/// @brief 2x2 matrix of real numbers
using Matrix2x2R = Matrix2x2<Real>;

/// @brief 3x3 matrix of real numbers
using Matrix3x3R = Matrix3x3<Real>;

/// @brief 3x3 matrix of real numbers
using Matrix4x4R = Matrix4x4<Real>;

/// @brief 2D vector of complex numbers
using Vector2C = Vector2<Complex>;

/// @brief 3D vector of complex numbers
using Vector3C = Vector3<Complex>;

/// @brief 4D vector of complex numbers
using Vector4C = Vector4<Complex>;

/// @brief 2x2 matrix of complex numbers
using Matrix2x2C = Matrix2x2<Complex>;

/// @brief 3x3 matrix of complex numbers
using Matrix3x3C = Matrix3x3<Complex>;

/// @brief 3x3 matrix of complex numbers
using Matrix4x4C = Matrix4x4<Complex>;

/// @}

}

#endif // _dc7a30fd_6048_4dfc_a551_efae434269f0
