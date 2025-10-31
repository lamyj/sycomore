#define BOOST_TEST_MODULE QuantityConstructors
#define BOOST_TEST_DYN_LINK

#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include "sycomore/Quantity.h"
#include "sycomore/QuantityArray.h"
#include "sycomore/QuantityTensor.h"
#include "sycomore/QuantityTensorFixed.h"

#include "utils.h"

using Types = boost::mpl::list<
    sycomore::Matrix2x2Q, sycomore::TensorQ<2>, sycomore::ArrayQ>;

using DynamicShapeTypes = boost::mpl::list<
    sycomore::TensorQ<2>, sycomore::ArrayQ>;

BOOST_AUTO_TEST_CASE_TEMPLATE(Shape, T, DynamicShapeTypes)
{
    T const q{typename T::Container::shape_type{2, 3}};
    BOOST_CHECK((q.shape() == typename T::Container::shape_type{2, 3}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ArrayOfScalars, T, Types)
{
    T const q{{1, 2}, {3, 4}};
    CHECK_QUANTITY(q, (T{{{1, 2}, {3, 4}}, sycomore::Dimensionless}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ArrayOfScalarsAndDimensions, T, Types)
{
    T const q{{{1, 2}, {3, 4}}, sycomore::Length};
    CHECK_QUANTITY(q, (T{{{1, 2}, {3, 4}}, sycomore::Length}));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(ArrayOfQuantity, T, Types)
{
    using sycomore::Quantity;
    using sycomore::Length;
    using sycomore::Time;
    using sycomore::Mass;
    using sycomore::ElectricCurrent;
    T const q{
        {Quantity(1, Length), Quantity(2, Length)},
        {Quantity(3, Length), Quantity(4, Length)}};
    CHECK_QUANTITY(q, (T{{{1, 2}, {3, 4}}, Length}));
    BOOST_CHECK_THROW(
        (T{
            {Quantity(1, Length), Quantity(1, Mass)},
            {Quantity(1, Time), Quantity(1, ElectricCurrent)}}),
        std::runtime_error);
}
