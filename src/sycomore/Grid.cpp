#include "Grid.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "sycomore/IndexGenerator.h"
#include "sycomore/magnetization.h"
#include "sycomore/sycomore.h"

namespace sycomore
{

Grid
::Grid()
: Grid(Index(), Shape())
{
    // Nothing else.
}

Grid
::Grid(Index const & origin, Shape const & shape)
: _origin(origin), _shape(shape), _stride(_compute_stride(shape))
{
    if(!this->_stride.empty())
    {
        this->_data.resize(this->_stride[this->_stride.size()-1]);
    }
}

Grid
::Grid(Index const & origin, Shape const & shape, value_type const & value)
: Grid(origin, shape)
{
    std::fill(this->_data.begin(), this->_data.end(), value);
}

Grid::value_type const *
Grid
::data() const
{
    return this->_data.data();
}

Grid::value_type *
Grid
::data()
{
    return this->_data.data();
}

Grid::value_type &
Grid
::operator[](Index const & index)
{
    Index translated_index(index.size());
    std::transform(
        index.begin(), index.end(), this->_origin.begin(),
        translated_index.begin(), std::minus<int>());
    auto const position = std::inner_product(
        translated_index.begin(), translated_index.end(), this->_stride.begin(),
        0);
    return this->_data[position];
}

Grid::value_type const &
Grid
::operator[](Index const & index) const
{
    Index translated_index(index.size());
    std::transform(
        index.begin(), index.end(), this->_origin.begin(),
        translated_index.begin(), std::minus<int>());
    auto const position = std::inner_product(
        translated_index.begin(), translated_index.end(), this->_stride.begin(),
        0);
    return this->_data[position];
}

size_t
Grid
::dimension() const
{
    return this->_shape.size();
}

Index const &
Grid
::origin() const
{
    return this->_origin;
}

Shape const &
Grid
::shape() const
{
    return this->_shape;
}

Stride const &
Grid
::stride() const
{
    return this->_stride;
}

void
Grid
::reshape(Index const & origin, Shape const & shape)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Grid new_grid(origin, shape);
    this->_reshape(new_grid);
}

void
Grid
::reshape(Index const & origin, Shape const & shape, value_type const & value)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Grid new_grid(origin, shape, value);
    this->_reshape(new_grid);
}

Grid::iterator
Grid
::begin()
{
    return this->_data.begin();
}

Grid::const_iterator
Grid
::begin() const
{
    return this->_data.begin();
}

Grid::const_iterator
Grid
::cbegin() const
{
    return this->_data.cbegin();
}

Grid::iterator
Grid
::end()
{
    return this->_data.end();
}

Grid::const_iterator
Grid
::end() const
{
    return this->_data.end();
}

Grid::const_iterator
Grid
::cend() const
{
    return this->_data.cend();
}

Stride
Grid::_compute_stride(Shape const & shape)
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

void
Grid
::_reshape(Grid & new_grid)
{
    Index max_origin(this->_origin.size());
    std::transform(
        this->_origin.begin(), this->_origin.end(), new_grid.origin().begin(),
        max_origin.begin(), [](int x, int y) { return std::max(x, y); });

    auto const old_last = IndexGenerator(this->origin(), this->shape()).last();
    auto const new_last = IndexGenerator(
        new_grid.origin(), new_grid.shape()).last();
    Index min_last(this->_origin.size());
    std::transform(
        old_last.begin(), old_last.end(), new_last.begin(),
        min_last.begin(), [](int x, int y) { return std::min(x, y); });

    // Check whether the two grids are disjoint (on at least one axis,
    // max_origin >= min_last)
    auto mismatch = std::mismatch(
        max_origin.begin(), max_origin.end(), min_last.begin(),
        std::less<int>());
    if(mismatch.first == max_origin.end())
    {
        Shape min_shape(min_last.size());
        std::transform(
            min_last.begin(), min_last.end(),
            max_origin.begin(), min_shape.begin(),
            std::minus<int>());

        // NOTE: we could memcopy line-by-line
        for(auto && index: IndexGenerator(max_origin, min_shape))
        {
            new_grid[index] = (*this)[index];
        }
    }
    // Otherwise the two grids are disjoint, nothing to copy.

    this->_origin = std::move(new_grid._origin);
    this->_shape = std::move(new_grid._shape);
    this->_stride = std::move(new_grid._stride);
    this->_data = std::move(new_grid._data);
}

}
