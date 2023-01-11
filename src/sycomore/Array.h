#ifndef _dc7a30fd_6048_4dfc_a551_efae434269f0
#define _dc7a30fd_6048_4dfc_a551_efae434269f0

#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

namespace sycomore
{

/// @brief Real numbers
using Real = double;

/// @brief Complex numbers
using Complex = std::complex<Real>;

/// @brief 3D Vector of arithmetic types
template <typename T>
using Vector3 = xt::xtensor_fixed<T, xt::xshape<3>>;

/// @brief 3x3 matrix of arithmetic types
template <typename T>
using Matrix3x3 = xt::xtensor_fixed<T, xt::xshape<3, 3>>;

/// @brief 4x4 matrix of arithmetic types
template <typename T>
using Matrix4x4 = xt::xtensor_fixed<T, xt::xshape<4, 4>>;

/// @brief Static-dimension array of real numbers
template<std::size_t N, xt::layout_type L=XTENSOR_DEFAULT_LAYOUT>
using TensorR = xt::xtensor<Real, N, L>;

/// @brief Dynamic-dimension array of real numbers
using ArrayR = xt::xarray<Real>;

/// @brief Dynamic-dimension array of complex numbers
using ArrayC = xt::xarray<Complex>;

/// @brief 3D vector of real numbers
using Vector3R = Vector3<Real>;

}

#endif // _dc7a30fd_6048_4dfc_a551_efae434269f0
