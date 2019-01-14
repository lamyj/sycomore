#include "Model.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/IndexGenerator.h"
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
    cl_device_type device_type)
: _species(species), _epsilon_squared(0)
{
    this->_initial_magnetization = as_complex_magnetization(magnetization);

    this->_time_intervals.resize(time_intervals.size());
    for(auto && item: time_intervals)
    {
        this->_time_intervals[this->_dimensions.size()] = item.second;
        this->_dimensions[item.first] = this->_dimensions.size();
    }

    Index const origin(time_intervals.size(), 0);
    Shape const shape (time_intervals.size(), 1);
    this->_bounding_box = {origin, shape};
    this->_m = Grid<ComplexMagnetization>(
        origin, shape, ComplexMagnetization(0,0,0));
    this->_m[Index(time_intervals.size(), 0)] = this->_initial_magnetization;

    Index F_origin(origin.size()+2, 0);
    std::copy(origin.begin(), origin.end(), F_origin.begin()+1);

    Shape F_shape(shape.size()+2, 0);
    std::copy(shape.begin(), shape.end(), F_shape.begin()+1);
    F_shape[0] = 3;
    F_shape[F_shape.size()-1] = this->_dimensions.size();

    this->_F = Grid<Real>(F_origin, F_shape, NAN);

    // Initialize OpenCL backend
    auto platforms = boost::compute::system::platforms();
    for(
        auto platform_it=platforms.begin();
        platform_it!=platforms.end() && !this->_device.id();
        ++platform_it)
    {
        auto && devices = platform_it->devices(device_type);
        for(auto && device: devices)
        {
            if(device.get_info<CL_DEVICE_DOUBLE_FP_CONFIG>())
            {
                this->_device = devices[0];
                break;
            }
        }
    }
    if(!this->_device.id())
    {
        throw std::runtime_error("No matching device");
    }

    this->_context = boost::compute::context(this->_device);
    this->_queue = boost::compute::command_queue(this->_context, this->_device);

    std::string source = BOOST_COMPUTE_STRINGIZE_SOURCE(
    inline double2 cmult(double2 c1, double2 c2)
    {
        return (double2)(c1.x*c2.x-c1.y*c2.y, c1.y*c2.x+c1.x*c2.y);
    }

    kernel void apply_pulse(
        global ComplexMagnetization * m_grid, global double2 const * R)
    {
        global ComplexMagnetization * m_ptr = m_grid+get_global_id(0);

        double2 p = cmult(R[0], m_ptr->p) + R[3]*m_ptr->z + cmult(R[6], m_ptr->m);
        double z = (cmult(R[1], m_ptr->p) + R[4]*m_ptr->z + cmult(R[7], m_ptr->m)).x;
        double2 m = cmult(R[2], m_ptr->p) + R[5]*m_ptr->z + cmult(R[8], m_ptr->m);

        m_ptr->p = p; m_ptr->z = z; m_ptr->m = m;
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

std::vector<TimeInterval> const &
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
    boost::compute::vector<ComplexMagnetization> m_d(
        this->_m.data(),
        this->_m.data()+this->_m.stride()[this->_m.dimension()],
        this->_queue);

    this->_apply_pulse.set_arg(0, m_d);
    this->_apply_pulse.set_arg(1, R_d);
    this->_queue.enqueue_1d_range_kernel(this->_apply_pulse, 0, m_d.size(), 0);

    boost::compute::copy(m_d.begin(), m_d.end(), this->_m.data(), this->_queue);
}

void
Model
::apply_time_interval(std::string const & name)
{
    auto && mu = this->_dimensions.at(name);
    auto && time_interval = this->_time_intervals[mu];

    // Update the bounding box, resize the grids if needed.
    --this->_bounding_box.first[mu];
    this->_bounding_box.second[mu] += 2;
    if(this->_bounding_box.first[mu] < this->_m.origin()[mu])
    {
        // Expand along mu
        auto origin = this->_m.origin();
        origin[mu] = std::min(-1, 2*origin[mu]-1);

        auto shape = this->_m.shape();
        shape[mu] = 2*shape[mu]+1;

        this->_m.reshape(origin, shape, ComplexMagnetization::zero);

        Index F_origin(this->_F.origin());
        std::copy(origin.begin(), origin.end(), F_origin.begin()+1);

        Shape F_shape(this->_F.shape());
        std::copy(shape.begin(), shape.end(), F_shape.begin()+1);
        F_shape[0] = 3;
        F_shape[F_shape.size()-1] = this->_dimensions.size();

        this->_F.reshape(F_origin, F_shape, NAN);
    }

    auto const E_1 = std::exp(-this->_species.R1 * time_interval.duration);
    auto const E_2 = std::exp(-this->_species.R2 * time_interval.duration);

    // Relaxation effects: update m- with the forward neighbor, m+ with the
    // backward neighbor and m0 with the current configuration.
    // Since the grids are kept symmetric, we can do it in a single, simple, scan.
    auto && m_stride = this->_m.stride()[mu];
    int const last = this->_m.origin()[mu]+this->_m.shape()[mu]-1;
    // Diffusion effects
    auto && p_mu = time_interval.gradient_moment;
    auto const p_mu_norm_third = dot(p_mu, p_mu)/3;
    auto const minus_D_tau = -this->_species.D*time_interval.duration;
    auto const do_diffusion =
        this->_species.D > 0 && p_mu != Array<Real>(p_mu.size(), 0);
    // WARNING: to use the offset, we need to iterate on the whole grid, not
    // on the bounding box only.
    std::pair<decltype(this->_m.data()), decltype(this->_m.data())> m_it;
    m_it.first = this->_m.data();
    m_it.second = m_it.first + this->_m.stride()[this->_m.dimension()] - 1;

    auto && F_stride = this->_F.stride()[1+mu];
    std::pair<decltype(this->_F.data()), decltype(this->_F.data())> F_it;
    F_it.first = this->_F.data() + mu*this->_F.stride()[this->_F.dimension()-1];
    F_it.second = F_it.first + this->_F.stride()[this->_F.dimension()-1] - 3;

    for(auto && index: IndexGenerator(this->_m.origin(), this->_m.shape()))
    {
        if(index[mu] != last)
        {
            Real F_minus = 1;
            Real F = 1;
            Real F_plus = 1;
            if(do_diffusion)
            {
                auto F_it_neighbor = F_it.first+F_stride;
                if(std::isnan(*F_it_neighbor))
                {
                    auto neighbor(index); ++neighbor[mu];
                    *(F_it_neighbor) = this->_compute_F(
                        neighbor, -1, p_mu, minus_D_tau, p_mu_norm_third);
                }
                F_minus = *(F_it_neighbor);

                if(std::isnan(*(F_it.first+1)))
                {
                    *(F_it.first+1) = this->_compute_F(
                        index, 0, p_mu, minus_D_tau, p_mu_norm_third);
                }
                F = *(F_it.first+1);

                F_it_neighbor = F_it.second-F_stride+2;
                if(std::isnan(*F_it_neighbor))
                {
                    auto neighbor(-index); --neighbor[mu];
                    *F_it_neighbor = this->_compute_F(
                        neighbor, +1, p_mu, minus_D_tau, p_mu_norm_third);
                }
                F_plus = *F_it_neighbor;
            }

            m_it.first->m = F_minus * E_2 * (m_it.first+m_stride)->m;
            m_it.first->z *= F * E_1;
            m_it.second->p = F_plus * E_2 * (m_it.second-m_stride)->p;
        }
        // Else no forward neighbor and thus no symmetric backward neighbor.

        ++m_it.first; --m_it.second;
        F_it.first += 3; F_it.second -= 3;
    }

    // Repolarization: second term of Equation 19
    auto const repolarization = this->_species.w * (1-E_1);
    this->_m[Index(this->_time_intervals.size(), 0)].p +=
        repolarization * this->_initial_magnetization.p;
    this->_m[Index(this->_time_intervals.size(), 0)].z +=
        repolarization * this->_initial_magnetization.z;
    this->_m[Index(this->_time_intervals.size(), 0)].m +=
        repolarization * this->_initial_magnetization.m;

    if(this->_epsilon_squared > 0)
    {
        this->_cleanup();
    }
}

Grid<ComplexMagnetization> const &
Model
::magnetization() const
{
    return this->_m;
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
        for(auto && m: this->_m)
        {
            update_isochromat(m);
        }
    }
    else
    {
        for(auto && configuration: configurations)
        {
            update_isochromat(this->_m[configuration]);
        }
    }

    return isochromat;
}

Real
Model
::_compute_F(
    Index const & n, int j,
    Array<Real> const & p_mu, Real minus_D_tau, Real p_mu_norm_third)
{
    Array<Real> p_n(3, 0);
    for(size_t d=0; d<this->_dimensions.size(); ++d)
    {
        auto && time_interval = this->_time_intervals[d];
        p_n += n[d] * time_interval.gradient_moment;
    }

    auto const F = std::exp(
        minus_D_tau
        *(dot(p_n, p_n) + j*dot(p_n, p_mu) + j*j*p_mu_norm_third));
    return F;
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

    Index const zero_i(this->_m.dimension(), 0);

    Index first(this->_m.dimension(), 0);
    Index last(this->_m.dimension(), 0);

    auto update_boundary = [&](Index const & index) {
        for(size_t d=0; d<this->_m.dimension(); ++d)
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

        auto && m = this->_m[index];
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
        for(size_t i=0; i<this->_m.dimension(); ++i)
        {
            int neighbors_count = 0;

            Index neighbor = index;

            neighbor[i] += 1;
            if(
                neighbor[i] < this->_m.origin()[i]+this->_m.shape()[i]
                && this->_m[neighbor] != ComplexMagnetization::zero)
            {
                ++neighbors_count;
            }

            neighbor[i] -= 2; // NOTE: go two steps back to skip the current index
            if(
                neighbor[i] >= this->_m.origin()[i]
                && this->_m[neighbor] != ComplexMagnetization::zero)
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

        this->_m[index] = ComplexMagnetization::zero;
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
