#ifndef _6a89e095_d9bf_4d67_b451_9e878d6797c0
#define _6a89e095_d9bf_4d67_b451_9e878d6797c0

#include "object_type_caster.h"

#include <algorithm>
#include <queue>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <xtensor/xfixed.hpp>

template<typename Container>
void assign_empty(std::vector<std::size_t> const & s, Container & c)
{
    c = Container(typename Container::shape_type{s.begin(), s.end()});
}

template<typename ET, typename FSH, xt::layout_type L>
void assign_empty(
    std::vector<std::size_t> const & s, xt::xtensor_fixed<ET, FSH, L> & c)
{
    c = xt::xtensor_fixed<ET, FSH, L>();
}

template<typename Container>
bool
object_type_caster<Container>
::load(pybind11::handle source, bool)
{
    using value_type = typename Container::value_type;
    using shape_type = typename Container::shape_type;
    
    if(pybind11::isinstance<pybind11::buffer>(source))
    {
        auto const info = source.cast<pybind11::buffer>().request();
        if(info.format != "O")
        {
            throw std::runtime_error("Must be an object buffer");
        }
        
        assign_empty(Self::_shape(info), value);
        this->_copy(info);
    }
    else if(pybind11::isinstance<pybind11::sequence>(source))
    {
        auto sequence = source.cast<pybind11::sequence>();
        assign_empty(Self::_shape(sequence), value);
        this->_copy(sequence);
    }
    else
    {
        throw std::runtime_error(
            "Cannot create array from "
            +pybind11::str(source).cast<std::string>());
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
    // for(std::size_t i=0; i<destination.ndim(); ++i)
    // {
    //     std::cout << strides[i] << " ";
    // }
    // for(std::size_t i=0; i<destination.ndim(); ++i)
    // {
    //     std::cout << destination.strides(i) << " ";
    // }
    // std::cout << std::endl;
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
    for(std::size_t i=0; i<info.size; ++i)
    {
        pybind11::handle handle(*source_pointer);
        PyObject * pointer = handle.ptr();
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
