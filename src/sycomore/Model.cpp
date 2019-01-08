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
    this->_initial_magnetization = as_complex_magnetization(magnetization);

    for(auto && item: time_intervals)
    {
        this->_time_intervals.insert(item);
        this->_dimensions[item.first] = this->_dimensions.size();
    }

    Index const origin(time_intervals.size(), 0);
    Shape const shape (time_intervals.size(), 1);
    this->_bounding_box = {origin, shape};
    this->_grid = Grid(origin, shape, ComplexMagnetization(0,0,0));
    this->_grid[Index(time_intervals.size(), 0)] = this->_initial_magnetization;

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

    std::string source = BOOST_COMPUTE_STRINGIZE_SOURCE(
    inline double2 mult_cc(double2 c1, double2 c2)
    {
        return (double2)(c1.x*c2.x-c1.y*c2.y, c1.y*c2.x+c1.x*c2.y);
    }

    inline double2 mult_cr(double2 c, double r)
    {
        return (double2)(c.x*r, c.y*r);
    }

    kernel void apply_pulse(
        global ComplexMagnetization * grid, global const double2 * R)
    {
        global ComplexMagnetization * m_ptr = grid+get_global_id(0);
        ComplexMagnetization m = *m_ptr;

        ComplexMagnetization result;

        result.p = mult_cc(R[0], m.p) + mult_cr(R[3], m.z) + mult_cc(R[6], m.m);
        result.z = (mult_cc(R[1], m.p) + mult_cr(R[4], m.z) + mult_cc(R[7], m.m)).x;
        result.m = mult_cc(R[2], m.p) + mult_cr(R[5], m.z) + mult_cc(R[8], m.m);

        *m_ptr = result;
    });
    source = boost::compute::type_definition<ComplexMagnetization>() + "\n" + source;

    auto program = boost::compute::program::build_with_source(
            source, this->_context);
    this->_apply_pulse = boost::compute::kernel(program, "apply_pulse");
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

    // WARNING: using only the bounding box takes more time that applying the
    // pulse on the whole array
    boost::compute::vector<ComplexMagnetization> grid_d(
        this->_grid.data(),
        this->_grid.data()+this->_grid.stride()[this->_grid.dimension()],
        this->_queue);

    this->_apply_pulse.set_arg(0, grid_d);
    this->_apply_pulse.set_arg(1, R_d);
    this->_queue.enqueue_1d_range_kernel(
        this->_apply_pulse, 0, grid_d.size(), 0);

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
    // Since the grid is kept symmetric, we can do it in a single, simple, scan.
    auto && stride = this->_grid.stride()[mu];
    int const last = this->_grid.origin()[mu]+this->_grid.shape()[mu]-1;
    auto * iterator = this->_grid.data();
    auto * reverse_iterator = iterator+this->_grid.stride()[this->_grid.dimension()]-1;
    // WARNING: to use the offset, we need to iterate on the whole grid, not
    // on the bounding box only.
    for(auto && index: IndexGenerator(this->_grid.origin(), this->_grid.shape()))
    {
        if(index[mu] == last)
        {
            // No forward neighbor and thus no symmetric backward neighbor.
            ++iterator;
            --reverse_iterator;
            continue;
        }

        auto && m_forward = *(iterator+stride);
        auto && m_backward = *(reverse_iterator-stride);

        iterator->m = E_2*m_forward.m;
        iterator->z *= E_1;
        reverse_iterator->p = E_2*m_backward.p;

        ++iterator;
        --reverse_iterator;
    }

    // Repolarization: second term of Equation 19
    auto const repolarization = this->_species.w * (1-E_1);
    this->_grid[Index(this->_time_intervals.size(), 0)].p +=
        repolarization * this->_initial_magnetization.p;
    this->_grid[Index(this->_time_intervals.size(), 0)].z +=
        repolarization * this->_initial_magnetization.z;
    this->_grid[Index(this->_time_intervals.size(), 0)].m +=
        repolarization * this->_initial_magnetization.m;

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
        auto const m_r = as_real_magnetization(m);
        isochromat.x += m_r.x;
        isochromat.y += m_r.y;
        isochromat.z += m_r.z;
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

    std::cout << "Bouding box shrunk to:";
    for(auto && x: this->_bounding_box.first) { std::cout << " " << x; };
    std::cout << ", ";
    for(auto && x: this->_bounding_box.second) { std::cout << " " << x; };
    std::cout << "\n";
}

}
