#ifndef _dc7a30fd_6048_4dfc_a551_efae434269f0
#define _dc7a30fd_6048_4dfc_a551_efae434269f0

#include <xtensor/xfixed.hpp>

#include <sycomore/hash.h>

namespace sycomore
{

template <typename T>
using Vector3 = xt::xtensor_fixed<T, xt::xshape<3>>;

template <typename T>
using Matrix3x3 = xt::xtensor_fixed<T, xt::xshape<3, 3>>;

}

#endif // _dc7a30fd_6048_4dfc_a551_efae434269f0
