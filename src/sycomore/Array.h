#ifndef _dc7a30fd_6048_4dfc_a551_efae434269f0
#define _dc7a30fd_6048_4dfc_a551_efae434269f0

#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

namespace sycomore
{

using Real = double;
using Complex = std::complex<Real>;

template <typename T>
using Vector3 = xt::xtensor_fixed<T, xt::xshape<3>>;

template <typename T>
using Matrix3x3 = xt::xtensor_fixed<T, xt::xshape<3, 3>>;

template <typename T>
using Matrix4x4 = xt::xtensor_fixed<T, xt::xshape<4, 4>>;

template<std::size_t N, xt::layout_type L=XTENSOR_DEFAULT_LAYOUT>
using TensorR = xt::xtensor<Real, N, L>;

using ArrayR = xt::xarray<Real>;

using ArrayC = xt::xarray<Complex>;

}

#endif // _dc7a30fd_6048_4dfc_a551_efae434269f0
