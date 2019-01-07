#include "Model.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/magnetization.h"
#include "sycomore/Pulse.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"

#include <iostream>

namespace sycomore
{

Model
::Model(
    Species const & species, Magnetization const & magnetization,
    std::vector<std::pair<std::string, TimeInterval>> const & time_intervals)
: _species(species), _dimensions(), _time_intervals(), _epsilon_squared(0)
{
    for(auto && item: time_intervals)
    {
        this->_time_intervals.insert(item);
        this->_dimensions[item.first] = this->_dimensions.size();
    }

    // Always keep a margin around the model to avoid checking the boundary
    // conditions
    Index const origin(time_intervals.size(), -2);
    Shape const shape (time_intervals.size(), 5);
    this->_bounding_box = {origin, shape};
    this->_grid = Grid(origin, shape, ComplexMagnetization(0,0,0));
    this->_grid[Index(time_intervals.size(), 0)] = as_complex_magnetization(magnetization);
}

std::map<std::string, size_t> const &
Model
::dimensions() const
{
    return this->_dimensions;
}

std::map<std::string, TimeInterval> const &
Model
::time_intervals() const
{
    return this->_time_intervals;
}

Real
Model
::epsilon() const
{
    return std::sqrt(this->_epsilon_squared);
}

void
Model
::epsilon(Real const & value)
{
    this->_epsilon_squared = std::pow(value, 2.);
}

void
Model
::apply_pulse(Pulse const & pulse)
{
    auto const actual_angle = this->B1 * pulse.angle;
    Pulse const actual_pulse{actual_angle, pulse.phase};
    auto const R = actual_pulse.rotation_matrix();

    for(auto & m: this->_grid)
    {
        ComplexMagnetization const result = {
            R[{0,0}]*m.plus + R[{0,1}]*m.zero + R[{0,2}]*m.minus,
            (R[{1,0}]*m.plus + R[{1,1}]*m.zero + R[{1,2}]*m.minus).real(),
            R[{2,0}]*m.plus + R[{2,1}]*m.zero + R[{2,2}]*m.minus
        };
        m = std::move(result);
    }
}

void
Model
::apply_time_interval(std::string const & name)
{
    auto const & time_interval = this->_time_intervals.at(name);

    Index const zero(this->_time_intervals.size(), 0);
    auto const mu = this->_dimensions.at(name);

    auto min_value = std::numeric_limits<Index::value_type>::max();
    auto max_value = std::numeric_limits<Index::value_type>::min();

    // Find configurations which participate in the update
    std::set<Index> indices;
    for(auto && index: IndexGenerator(this->_grid.origin(), this->_grid.shape()))
    {
        auto && m = this->_grid[index];
        if(std::norm(m.plus) > 0)
        {
            auto n = index;
            ++n[mu];
            indices.insert(n);

            min_value = std::min(min_value, n[mu]);
            max_value = std::max(max_value, n[mu]);
        }
        if(m.zero > 0)
        {
            indices.insert(index);

            min_value = std::min(min_value, index[mu]);
            max_value = std::max(max_value, index[mu]);
        }
        if(std::norm(m.minus) > 0)
        {
            auto n = index;
            --n[mu];
            indices.insert(n);

            min_value = std::min(min_value, n[mu]);
            max_value = std::max(max_value, n[mu]);
        }
    }

    this->_bounding_box.first[mu] = std::min(
        this->_bounding_box.first[mu], min_value);
    this->_bounding_box.second[mu] = std::max(
        this->_bounding_box.second[mu],
        static_cast<Shape::value_type>(max_value-min_value+1));

    auto boundaries = std::minmax_element(
        indices.begin(), indices.end(),
        [&](Index const & x, Index const & y) { return x[mu] < y[mu]; });
    auto boundary = std::max(
        std::abs((*(boundaries.first))[mu]),
        std::abs((*(boundaries.second))[mu]));
    auto origin = this->_grid.origin();
    auto shape = this->_grid.shape();
    auto const radius = -origin[mu];

    // Always keep a margin around the model to avoid checking the boundary
    // conditions
    if(boundary > radius-2)
    {
        // Expand along mu
        origin[mu] = std::min(-1, 2*origin[mu]-1);
        shape[mu] = 2*shape[mu]+1;
    }

    Grid new_grid(origin, shape, ComplexMagnetization(0,0,0));

    auto const E_1 = std::exp(-this->_species.R1 * time_interval.duration);
    auto const E_2 = std::exp(-this->_species.R2 * time_interval.duration);

    static ComplexMagnetization const m0(0,0,0);

    for(auto && index: indices)
    {
        auto n = index;
        // No need to check the boundary conditions: we always have a margin
        // around the model.
        --n[mu]; auto && m_minus = this->_grid[n];
        ++n[mu]; auto && m = this->_grid[n];
        ++n[mu]; auto && m_plus = this->_grid[n];
        new_grid[index] = {E_2*m_minus.plus, E_1*m.zero, E_2*m_plus.minus};
    }

    // Repolarization: second term of Equation 19
    // WARNING: this assumes m_eq = [0,0,1] (as in CoMoTk)
    new_grid[zero].zero += this->_species.w * (1-E_1);

    if(this->_epsilon_squared > 0)
    {
        this->_cleanup(new_grid, this->_bounding_box);
    }
    this->_grid = std::move(new_grid);
}

Grid const &
Model
::grid() const
{
    return this->_grid;
}

Magnetization
Model
::isochromat(std::vector<Index> const & configurations) const
{
    Magnetization isochromat{0, 0, 0};

    auto const update_isochromat = [&](ComplexMagnetization const & m) {
        isochromat.x += m.plus.real();
        isochromat.y += m.plus.imag();
        isochromat.z += m.zero;
    };

    if(configurations.empty())
    {
        for(auto && m: this->_grid)
        {
            update_isochromat(m);
        }
    }
    else
    {
        for(auto && configuration: configurations)
        {
            update_isochromat(this->_grid[configuration]);
        }
    }

    isochromat.x *= std::sqrt(Real(2));
    isochromat.y *= std::sqrt(Real(2));

    return isochromat;
}

void
Model
::_cleanup(Grid & model, std::pair<Index, Shape> const & bounding_box)
{
    /*
     * A bounding box of the occupied configurations can be constructed
     * iteratively: at the starting configuration, the bounding box is [0]. When
     * applying the mu^th time interval, we compute the set of new
     * configurations and thus can update the bounding box on the mu-th axis.
     * NOTE: pulses do not change the set of occupied configurations
     *
     * We assume that the current state is convex: it is true at the start, and
     * true after applying a time interval since this action will "spread" the
     * magnetization but not create new empty states.
     *
     * Generate the set of indices spanning this bounding box and sort them
     * according to their distance (TODO: which one?) to the origin.
     *
     * For each configuration
     *   If it is 0
     *     Continue
     *   If its population is large enough
     *     Continue
     *   For each dimension
     *     If the configuration has two non-zero neighbors along this axis
     *       Continue to next configuration (its removal would create a concavity)
     */

    IndexGenerator const generator(bounding_box.first, bounding_box.second);
    std::vector<Index> indices(generator.begin(), generator.end());

    Index const zero_i(model.dimension(), 0);
    ComplexMagnetization const zero_m(0, 0, 0);

    // Careful with OpenCL here: could removing two points in parallel lead to
    // a concavity?
    for(auto && index: indices)
    {
        if(index == zero_i)
        {
            continue;
        }

        auto && m = model[index];
        auto const magnitude =
            m.plus * std::conj(m.plus)
            + m.zero*m.zero
            + m.minus * std::conj(m.minus);
        if(magnitude.real() >= this->_epsilon_squared)
        {
            continue;
        }

        bool would_create_concavity = false;
        for(size_t i=0; i<model.dimension(); ++i)
        {
            int neighbors_count = 0;

            Index neighbor = index;

            neighbor[i] += 1;
            if(neighbor[i] < model.origin()[i]+model.shape()[i] && model[neighbor] != zero_m)
            {
                ++neighbors_count;
            }

            neighbor[i] -= 2; // NOTE: go two steps back to skip the current index
            if(neighbor[i] >= model.origin()[i] && model[neighbor] != zero_m)
            {
                ++neighbors_count;
            }

            if(neighbors_count == 2)
            {
                would_create_concavity = true;
                break;
            }
        }
        if(would_create_concavity)
        {
            continue;
        }

        model[index] = zero_m;
    }

}

}
