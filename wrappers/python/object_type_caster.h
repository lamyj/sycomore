#ifndef _50132908_42b6_425f_8bbe_21054c5c49f9
#define _50132908_42b6_425f_8bbe_21054c5c49f9

#include <string>
#include <type_traits>
#include <vector>

#include <pybind11/pybind11.h>

#if PYBIND11_VERSION_HEX <= 0x02090100
namespace
{
using std::enable_if_t;
using std::remove_cv_t;
using namespace pybind11;
}
#endif

template<typename Container>
struct object_type_caster
{
public:
    using Self = object_type_caster<Container>;
    
    PYBIND11_TYPE_CASTER(
        Container, pybind11::detail::_("numpy.ndarray[")
            + pybind11::detail::type_caster<typename Container::value_type>::name
            + pybind11::detail::_("]"));
    
    bool load(pybind11::handle source, bool);
    
    static pybind11::handle cast(
        Container const & source,
        pybind11::return_value_policy, pybind11::handle);
    
private:
    static std::vector<std::size_t> _shape(pybind11::buffer_info const & info);
    static std::vector<std::size_t> _shape(pybind11::sequence const & sequence);
    
    void _copy(pybind11::buffer_info const & info);
    void _copy(pybind11::sequence sequence);
};

#include "object_type_caster.txx"

#endif // _50132908_42b6_425f_8bbe_21054c5c49f9
