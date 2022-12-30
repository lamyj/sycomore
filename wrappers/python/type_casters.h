#ifndef _20f638ed_34c7_4438_824b_0dba29bf2694
#define _20f638ed_34c7_4438_824b_0dba29bf2694

#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

#include <xtensor-python/pyarray.hpp>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/Quantity.h"

#include "object_type_caster.h"

namespace pybind11
{

namespace detail
{

template<xt::layout_type L>
struct type_caster<xt::xarray<sycomore::Quantity, L>>:
    public object_type_caster<xt::xarray<sycomore::Quantity, L>>
{
};

template<std::size_t N, xt::layout_type L>
struct type_caster<xt::xtensor<sycomore::Quantity, N, L>>:
    public object_type_caster<xt::xtensor<sycomore::Quantity, N, L>>
{
};

template<typename FSH, xt::layout_type L>
struct type_caster<xt::xtensor_fixed<sycomore::Quantity, FSH, L>>:
    public object_type_caster<xt::xtensor_fixed<sycomore::Quantity, FSH, L>>
{
};

}
}

#endif // _20f638ed_34c7_4438_824b_0dba29bf2694
