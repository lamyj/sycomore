#ifndef _e55feb28_6959_4759_ac43_7eaddb0fb7a3
#define _e55feb28_6959_4759_ac43_7eaddb0fb7a3

#include <boost/test/tools/assertion_result.hpp>

#include <xtensor/xmath.hpp>

template<typename T>
auto difference(T const & left, T const & right, double epsilon=1e-12)
{
    return xt::amax(xt::abs(left - right))();
}

inline double difference(double left, double right, double epsilon=1e-12)
{
    return std::abs(left-right);
}

template<typename T>
boost::test_tools::predicate_result
is_close(T const & left, T const & right, double epsilon=1e-12)
{
    boost::test_tools::predicate_result result(true);
    auto const delta = difference(left, right, epsilon);
    if(delta >= epsilon)
    {
        result = false;
        result.message()
            << "Difference exceeds tolerance "
            << "[" << delta << " > " << epsilon << "]";
    }
    return result;
}

#define CHECK_IDENTITY(l, r) \
    BOOST_CHECK(typeid(decltype(l)) == typeid(decltype(r))); \
    BOOST_CHECK(&l == &r);

#define CHECK_QUANTITY(l, r) \
    BOOST_CHECK(is_close((l).magnitude, (r).magnitude)); \
    BOOST_CHECK((l).dimensions == (r).dimensions);

#define CHECK_TYPE_AND_QUANTITY(x, T, z) \
    BOOST_CHECK(typeid(decltype(x)) == typeid(T)); \
    CHECK_QUANTITY(x, z);

#define CHECK_IDENTITY_AND_QUANTITY(x, y, z) \
    CHECK_IDENTITY(x, y); \
    CHECK_QUANTITY(x, z);

#endif // _e55feb28_6959_4759_ac43_7eaddb0fb7a3
