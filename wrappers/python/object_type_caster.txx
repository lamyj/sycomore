#ifndef _6a89e095_d9bf_4d67_b451_9e878d6797c0
#define _6a89e095_d9bf_4d67_b451_9e878d6797c0

#include "object_type_caster.h"

#include <algorithm>
#include <queue>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

// NOTE: the initialization of the shape type differs across containers. We need
// to have specialization for each container type.

template<typename ET, xt::layout_type L>
bool assign_empty(std::vector<std::size_t> const & s, xt::xarray<ET, L> & c)
{
    typename xt::xarray<ET, L>::shape_type const shape(s.begin(), s.end());
    c = xt::xarray<ET, L>(shape);
    return true;
}

template<typename ET, std::size_t N, xt::layout_type L>
bool assign_empty(std::vector<std::size_t> const & s, xt::xtensor<ET, N, L> & c)
{
    if(s.size() != N)
    {
        return false;
    }
    
    typename xt::xtensor<ET, N, L>::shape_type shape;
    std::copy(s.begin(), s.end(), shape.begin());
    c = xt::xtensor<ET, N, L>(shape);
    return true;
}


template<typename ET, typename FSH, xt::layout_type L>
bool assign_empty(
    std::vector<std::size_t> const & s, xt::xtensor_fixed<ET, FSH, L> & c)
{
    if(!xt::same_shape(s, FSH()))
    {
        return false;
    }
    c = xt::xtensor_fixed<ET, FSH, L>();
    return true;
}

template<typename Container>
bool
object_type_caster<Container>
::load(pybind11::handle source, bool)
{
    if(pybind11::isinstance<pybind11::buffer>(source))
    {
        auto const info = source.cast<pybind11::buffer>().request();
        if(info.format != "O")
        {
            return false;
        }
        
        std::vector<std::size_t> shape;
        try
        {
            shape = Self::_shape(info);
        }
        catch(std::runtime_error const &)
        {
            return false;
        }
        if(!assign_empty(shape, value))
        {
            return false;
        }
        this->_copy(info);
    }
    else if(pybind11::isinstance<pybind11::sequence>(source))
    {
        auto sequence = source.cast<pybind11::sequence>();
        
        std::vector<std::size_t> shape;
        try
        {
            shape = Self::_shape(sequence);
        }
        catch(std::runtime_error const &)
        {
            return false;
        }
        
        if(!assign_empty(shape, value))
        {
            return false;
        }
        this->_copy(sequence);
    }
    else
    {
        return false;
    }
    
    return (PyErr_Occurred() == NULL);
}

template<typename Container>
pybind11::handle
object_type_caster<Container>
::cast(
    Container const & source,
    pybind11::return_value_policy, pybind11::handle)
{
    auto source_ptr = source.data();
    
    // WARNING: the strides must be in size of destination type, i.e. PyObject*
    std::vector<std::size_t> strides(
        source.strides().begin(), source.strides().end());
    std::for_each(
        strides.begin(), strides.end(),
        [](auto & x){ x *= sizeof(PyObject*); });
    pybind11::array destination(pybind11::dtype("O"), source.shape(), strides);
    
    auto destination_ptr = reinterpret_cast<PyObject**>(
        destination.mutable_data());
    
    for(std::size_t i=0; i<source.size(); ++i)
    {
        *destination_ptr = pybind11::cast(*source_ptr).release().ptr();
        ++source_ptr;
        ++destination_ptr;
    }
    return destination.release();
}

template<typename Container>
std::vector<std::size_t>
object_type_caster<Container>
::_shape(pybind11::buffer_info const & info)
{
    return {info.shape.begin(), info.shape.end()};
}

template<typename Container>
std::vector<std::size_t>
object_type_caster<Container>
::_shape(pybind11::sequence const & sequence)
{
    // Get the shape of each subsequence item
    std::vector<std::vector<std::size_t>> shapes;
    shapes.reserve(sequence.size());
    for(auto && item: sequence)
    {
        if(pybind11::isinstance<pybind11::sequence>(item))
        {
            shapes.push_back(Self::_shape(item));
        }
    }
    
    if(shapes.size() == 0)
    {
        // Sequence of scalars
        return {sequence.size()};
    }
    else if(shapes.size() != sequence.size())
    {
        // Mix of scalars and subsequences
        throw std::runtime_error("Invalid nesting (mixed sequence)");
    }
    else
    {
        // Nested sequence
        auto const subshape = shapes[0];
        
        // Check whether all shapes agree
        auto const correct_nesting = std::all_of(
            shapes.begin(), shapes.end(), [&](auto && x) {return x==subshape;});
        if(!correct_nesting)
        {
            throw std::runtime_error("Invalid nesting (mixed shapes)");
        }
        
        // Build final shape
        std::vector<std::size_t> shape(1+subshape.size());
        shape[0] = sequence.size();
        std::copy(subshape.begin(), subshape.end(), shape.begin()+1);
        return shape;
    }
}

template<typename Container>
void
object_type_caster<Container>
::_copy(pybind11::buffer_info const & info)
{
    auto source_pointer = reinterpret_cast<PyObject**>(info.ptr);
    auto destination_pointer = this->value.data();
    for(std::ptrdiff_t i=0; i<info.size; ++i)
    {
        pybind11::handle handle(*source_pointer);
        *destination_pointer = handle.cast<typename Container::value_type>();
        ++source_pointer;
        ++destination_pointer;
    }
}

template<typename Container>
void
object_type_caster<Container>
::_copy(pybind11::sequence sequence)
{
    std::queue<pybind11::sequence> queue({sequence});
    for(std::size_t i=0; i<this->value.dimension()-1; ++i)
    {
        for(std::size_t n=0, end=queue.size(); n<end; ++n)
        {
            auto && item = queue.front();
            for(auto && child: item)
            {
                queue.push(child);
            }
            queue.pop();
        }
    }
    
    auto destination = this->value.begin();
    while(!queue.empty())
    {
        auto && item = queue.front();
        destination = std::transform(
            item.begin(), item.end(), destination,
            [](auto && x) {
                return x.template cast<typename Container::value_type>(); });
        queue.pop();
    }
}

#endif // _6a89e095_d9bf_4d67_b451_9e878d6797c0
