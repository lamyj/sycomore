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

BOOST_COMPUTE_ADAPT_STRUCT(
    sycomore::ComplexMagnetization, ComplexMagnetization, (p, z, m))

namespace sycomore
{

Model
::Model(
    Species const & species, Magnetization const & magnetization,
    std::vector<std::pair<std::string, TimeInterval>> const & time_intervals,
    cl_device_type device_filter)
: _species(species), _epsilon_squared(0)
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

    // Initialize OpenCL backend
    for(auto && platform: boost::compute::system::platforms())
    {
        auto && devices = platform.devices(
            device_filter & CL_DEVICE_DOUBLE_FP_CONFIG);
        if(!devices.empty())
        {
            this->_device = devices[0];
            break;
        }
    }
    if(!this->_device.id())
    {
        throw std::runtime_error("No matching device");
    }

    this->_context = boost::compute::context(this->_device);
    this->_queue = boost::compute::command_queue(this->_context, this->_device);
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

    boost::compute::vector<Complex> R_d(
        R.data(), R.data()+R.stride()[R.dimension()], this->_queue);
    BOOST_COMPUTE_CLOSURE(
        ComplexMagnetization, apply_pulse, (ComplexMagnetization m), (R_d),
        {
            ComplexMagnetization result;
            result.p.x =
                R_d[0].x*m.p.x - R_d[0].y*m.p.y
                + R_d[3].x*m.z
                + R_d[6].x*m.m.x - R_d[6].y*m.m.y;
            result.p.y =
                R_d[0].x*m.p.y + R_d[0].y*m.p.x
                + R_d[3].y*m.z
                + R_d[6].x*m.m.y + R_d[6].y*m.m.x;

            result.z =
                R_d[1].x*m.p.x - R_d[1].y*m.p.y
                + R_d[4].x*m.z
                + R_d[7].x*m.m.x - R_d[7].y*m.m.y;

            result.m.x =
                R_d[0].x*m.p.x - R_d[0].y*m.p.y
                + R_d[3].x*m.z
                + R_d[6].x*m.m.x - R_d[6].y*m.m.y;
            result.m.y =
                R_d[2].x*m.p.y + R_d[2].y*m.p.x
                + R_d[5].y*m.z
                + R_d[8].x*m.m.y + R_d[8].y*m.m.x;
            return result;
        });
    boost::compute::vector<ComplexMagnetization> grid_d(
        this->_grid.data(), this->_grid.data()+this->_grid.stride()[this->_grid.dimension()],
        this->_queue);
    boost::compute::transform(
        grid_d.begin(), grid_d.end(), grid_d.begin(), apply_pulse,
        this->_queue);
    boost::compute::copy(
        grid_d.begin(), grid_d.end(), this->_grid.data(), this->_queue);
}

void
Model
::apply_time_interval(std::string const & name)
{
    auto && time_interval = this->_time_intervals.at(name);
    auto && mu = this->_dimensions.at(name);

    // Update the bounding box, resize the grid if needed.
    --this->_bounding_box.first[mu];
    this->_bounding_box.second[mu] += 2;
    if(this->_bounding_box.first[mu] < this->_grid.origin()[mu])
    {
        // Expand along mu
        auto origin = this->_grid.origin();
        origin[mu] = std::min(-1, 2*origin[mu]-1);

        auto shape = this->_grid.shape();
        shape[mu] = 2*shape[mu]+1;

        this->_grid.reshape(origin, shape, ComplexMagnetization::zero);
    }

    auto const E_1 = std::exp(-this->_species.R1 * time_interval.duration);
    auto const E_2 = std::exp(-this->_species.R2 * time_interval.duration);

    // Relaxation effects: update m- with the forward neighbor, m+ with the
    // backward neighbor and m0 with the current configuration.
    size_t offset = 0;
    auto * data = this->_grid.data();
    auto const stride = this->_grid.stride()[mu];
    // WARNING: to use the offset, we need to iterate on the whole grid, not
    // on the bounding box only.
    for(auto && index: IndexGenerator(this->_grid.origin(), this->_grid.shape()))
    {
        if(index[mu] == this->_grid.origin()[mu]+this->_grid.shape()[mu]-1)
        {
            // No forward neighbor and no symmetric backward neighbor.
            ++offset;
            continue;
        }

        auto && m_forward = data[offset+stride];

        auto const symmetric_offset = this->_grid.stride()[this->_grid.dimension()]-offset-1;
        auto && m_backward = data[symmetric_offset-stride];

        data[offset].m = E_2*m_forward.m;
        data[offset].z *= E_1;
        data[symmetric_offset].p = E_2*m_backward.p;

        ++offset;
    }

    // Repolarization: second term of Equation 19
    // WARNING: this assumes m_eq = [0,0,1] (as in CoMoTk)
    new_grid[zero].zero += this->_species.w * (1-E_1);

    if(this->_epsilon_squared > 0)
    {
        this->_cleanup();
    }
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
        isochromat.x += m.p.real();
        isochromat.y += m.p.imag();
        isochromat.z += m.z;
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
::_cleanup()
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

    Index const zero_i(this->_grid.dimension(), 0);

    Index first(this->_grid.dimension(), 0);
    Index last(this->_grid.dimension(), 0);

    auto update_boundary = [&](Index const & index) {
        for(size_t d=0; d<this->_grid.dimension(); ++d)
        {
            first[d] = std::min(first[d], index[d]);
            last[d] = std::max(first[d], index[d]);
        }
    };

    for(auto && index: IndexGenerator(this->_bounding_box.first, this->_bounding_box.second))
    {
        if(index == zero_i)
        {
            update_boundary(index);
            continue;
        }

        auto && m = this->_grid[index];
        auto const magnitude =
            m.p * std::conj(m.p)
            + m.z*m.z
            + m.m * std::conj(m.m);
        if(magnitude.real() >= this->_epsilon_squared)
        {
            update_boundary(index);
            continue;
        }

        bool would_create_concavity = false;
        for(size_t i=0; i<this->_grid.dimension(); ++i)
        {
            int neighbors_count = 0;

            Index neighbor = index;

            neighbor[i] += 1;
            if(
                neighbor[i] < this->_grid.origin()[i]+this->_grid.shape()[i]
                && this->_grid[neighbor] != ComplexMagnetization::zero)
            {
                ++neighbors_count;
            }

            neighbor[i] -= 2; // NOTE: go two steps back to skip the current index
            if(
                neighbor[i] >= this->_grid.origin()[i]
                && this->_grid[neighbor] != ComplexMagnetization::zero)
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
            update_boundary(index);
            continue;
        }

        this->_grid[index] = ComplexMagnetization::zero;
    }

    Shape shape(first.size());
    std::transform(
        first.begin(), first.end(), last.begin(), shape.begin(),
        [](int f, int l) { return l-f+1;});
    this->_bounding_box = {first, shape};
}

}
