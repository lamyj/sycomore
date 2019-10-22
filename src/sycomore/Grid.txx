#include "Grid.h"

#include <algorithm>
#include <utility>
#include <vector>

#include "sycomore/GridScanner.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

template<typename T>
Grid<T>
::Grid()
: Grid(Index(), Shape())
{
    // Nothing else.
}

template<typename T>
Grid<T>
::Grid(Index const & origin, Shape const & shape)
: _origin(origin), _shape(shape), _stride(_compute_stride(shape))
{
    if(!this->_stride.empty())
    {
        this->_data.resize(this->_stride[this->_stride.size()-1]);
    }
}

template<typename T>
Grid<T>
::Grid(Index const & origin, Shape const & shape, value_type const & value)
: Grid(origin, shape)
{
    std::fill(this->_data.begin(), this->_data.end(), value);
}

template<typename T>
typename Grid<T>::value_type const *
Grid<T>
::data() const
{
    return this->_data.data();
}

template<typename T>
typename Grid<T>::value_type *
Grid<T>
::data()
{
    return this->_data.data();
}

template<typename T>
typename Grid<T>::value_type &
Grid<T>
::operator[](std::size_t offset)
{
    return this->_data[offset];
}

template<typename T>
typename Grid<T>::value_type const &
Grid<T>
::operator[](std::size_t offset) const
{
    return this->_data[offset];
}

template<typename T>
typename Grid<T>::value_type &
Grid<T>
::operator[](Index const & index)
{
    auto const offset = dot(index-this->_origin, this->_stride);
    return this->operator[](offset);
}

template<typename T>
typename Grid<T>::value_type const &
Grid<T>
::operator[](Index const & index) const
{
    auto const offset = dot(index-this->_origin, this->_stride);
    return this->operator[](offset);
}

template<typename T>
std::size_t
Grid<T>
::dimension() const
{
    return this->_shape.size();
}

template<typename T>
Index const &
Grid<T>
::origin() const
{
    return this->_origin;
}

template<typename T>
Shape const &
Grid<T>
::shape() const
{
    return this->_shape;
}

template<typename T>
Stride const &
Grid<T>
::stride() const
{
    return this->_stride;
}

template<typename T>
void
Grid<T>
::reshape(Index const & origin, Shape const & shape)
{
    Grid new_grid(origin, shape);
    this->_reshape(new_grid);
}

template<typename T>
void
Grid<T>
::reshape(Index const & origin, Shape const & shape, value_type const & value)
{
    Grid new_grid(origin, shape, value);
    this->_reshape(new_grid);
}

template<typename T>
typename Grid<T>::iterator
Grid<T>
::begin()
{
    return this->_data.begin();
}

template<typename T>
typename Grid<T>::const_iterator
Grid<T>
::begin() const
{
    return this->_data.begin();
}

template<typename T>
typename Grid<T>::const_iterator
Grid<T>
::cbegin() const
{
    return this->_data.cbegin();
}

template<typename T>
typename Grid<T>::iterator
Grid<T>
::end()
{
    return this->_data.end();
}

template<typename T>
typename Grid<T>::const_iterator
Grid<T>
::end() const
{
    return this->_data.end();
}

template<typename T>
typename Grid<T>::const_iterator
Grid<T>
::cend() const
{
    return this->_data.cend();
}

template<typename T>
Stride
Grid<T>
::_compute_stride(Shape const & shape)
{
    if(shape.size() == 0)
    {
        return Stride();
    }

    Stride stride(shape.size()+1);
    stride[0] = 1;
    for(unsigned int i=0; i< shape.size(); ++i)
    {
        stride[i+1] = shape[i]*stride[i];
    }

    return stride;
}

template<typename T>
void
Grid<T>
::_reshape(Grid & new_grid)
{
    if(this->dimension() != new_grid.dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Index const max_origin = maximum(this->_origin, new_grid.origin());

    auto const old_end = GridScanner(
        this->origin(), this->shape()).region_end();
    auto const new_end = GridScanner(
        new_grid.origin(), new_grid.shape()).region_end();
    Index min_end = minimum(old_end, new_end);

    // Check whether the two grids are disjoint (on at least one axis,
    // max_origin >= min_last)
    auto mismatch = std::mismatch(
        max_origin.begin(), max_origin.end(), min_end.begin(),
        std::less<int>());
    if(mismatch.first == max_origin.end())
    {
        Shape min_shape(min_end.size());
        std::transform(
            min_end.begin(), min_end.end(),
            max_origin.begin(), min_shape.begin(),
            std::minus<int>());

        GridScanner const source_scanner(
            this->_origin, this->_shape, max_origin, min_shape);
        auto source_it = source_scanner.begin();
        auto const source_end = source_scanner.end();
        GridScanner const destination_scanner(
            new_grid._origin, new_grid._shape, max_origin, min_shape);
        auto destination_it = destination_scanner.begin();
        auto const destination_end = destination_scanner.end();
        while(source_it != source_end && destination_it != destination_end)
        {
            new_grid._data[destination_it->second] = this->_data[source_it->second];
            ++source_it;
            ++destination_it;
        }
    }
    // Otherwise the two grids are disjoint, nothing to copy.

    this->_origin = std::move(new_grid._origin);
    this->_shape = std::move(new_grid._shape);
    this->_stride = std::move(new_grid._stride);
    this->_data = std::move(new_grid._data);
}

}
