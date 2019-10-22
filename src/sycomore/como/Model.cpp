#include "Model.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/GridScanner.h"
#include "sycomore/magnetization.h"
#include "sycomore/Pulse.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"

namespace sycomore
{

namespace como
{

Model
::Model(
    Species const & species, Magnetization const & magnetization,
    std::vector<std::pair<std::string, TimeInterval>> const & time_intervals)
: _species(species), _epsilon_squared(0)
{
    this->_initial_magnetization = as_complex_magnetization(magnetization);

    this->_time_intervals.resize(time_intervals.size());
    for(auto && item: time_intervals)
    {
        this->_time_intervals[this->_dimensions.size()] = item.second;
        this->_dimensions.emplace(item.first, this->_dimensions.size());
    }

    // Grid of magnetization (i.e. configuration vectors)
    Index const origin(time_intervals.size(), 0);
    Shape const shape (time_intervals.size(), 1);

    this->_m = Grid<ComplexMagnetization>(
        origin, shape, ComplexMagnetization(0,0,0));
    this->_m[Index(time_intervals.size(), 0)] = this->_initial_magnetization;

    this->_bounding_box = {origin, shape};

    // Grid of time scales
    this->_tau = Grid<Real>(this->_m.origin(), this->_m.shape(), NAN);
    this->_tau[Index(origin.size(), 0)] = 0;

    // Grid of gradient moments
    Index p_origin(origin.size()+1, 0);
    std::copy(origin.begin(), origin.end(), p_origin.begin()+1);

    Shape p_shape(shape.size()+1, 0);
    std::copy(shape.begin(), shape.end(), p_shape.begin()+1);
    p_shape[0] = 3;

    this->_p = Grid<Real>(p_origin, p_shape, NAN);
    Array<Real> p(this->_p.data()+dot(-this->_p.origin(), this->_p.stride()), 3);
    std::fill(p.begin(), p.end(), 0);

    // Grid of diffusion damping factors
    Index F_origin(origin.size()+2, 0);
    std::copy(origin.begin(), origin.end(), F_origin.begin()+1);

    Shape F_shape(shape.size()+2, 0);
    std::copy(shape.begin(), shape.end(), F_shape.begin()+1);
    F_shape[0] = 3;
    F_shape[F_shape.size()-1] = this->_dimensions.size();

    this->_F = Grid<Real>(F_origin, F_shape, NAN);
}

std::map<std::string, std::size_t> const &
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
::get_epsilon() const
{
    return std::sqrt(this->_epsilon_squared);
}

void
Model
::set_epsilon(Real const & value)
{
    this->_epsilon_squared = std::pow(value, 2.);
}

void
Model
::apply_pulse(Pulse const & pulse)
{
    auto const start = std::chrono::high_resolution_clock::now();

    auto const R_m = pulse.rotation_matrix();
    auto && R = R_m.data();

    #pragma omp parallel for
#ifdef _WIN32
    // WARNING: only signed integer types in OpenMP loops on Windows
    for(int index=0; index<this->_m.stride()[this->_m.dimension()]; ++index)
    {
        auto it = this->_m.data()+index;
#else
    for(
        auto it=this->_m.data();
        it<this->_m.data()+this->_m.stride()[this->_m.dimension()]; ++it)
    {
#endif
        auto && m = *it;

        ComplexMagnetization const result{
            R[0]*m.p + R[3]*m.z + R[6]*m.m,
            (R[1]*m.p + R[4]*m.z + R[7]*m.m).real(),
            R[2]*m.p + R[5]*m.z + R[8]*m.m};

        *it = std::move(result);
    }

    this->_timers["pulse"] +=
        std::chrono::duration<double, std::ratio<1>>(
            std::chrono::high_resolution_clock::now()-start).count();
}

void
Model
::apply_pulse(HardPulseApproximation const & pulse)
{
    if(pulse.get_pulses().empty())
    {
        return;
    }

    this->apply_pulse(pulse.get_pulses()[0]);
    for(std::size_t i=1; i!=pulse.get_pulses().size(); ++i)
    {
        this->apply_time_interval(pulse.get_name());
        this->apply_pulse(pulse.get_pulses()[i]);
    }
}

void
Model
::apply_time_interval(std::string const & name)
{
    using namespace units;

    auto const start = std::chrono::high_resolution_clock::now();

    auto && mu = this->_dimensions.at(name);
    auto && time_interval = this->_time_intervals[mu];

    // Update the bounding box, resize the grids if needed.
    --this->_bounding_box.first[mu];
    this->_bounding_box.second[mu] += 2;
    if(this->_bounding_box.first[mu] < this->_m.origin()[mu])
    {
        // Expand m along mu
        auto origin = this->_m.origin();
        origin[mu] = std::min(-1, 2*origin[mu]-1);

        auto shape = this->_m.shape();
        shape[mu] = 2*shape[mu]+1;

        this->_m.reshape(origin, shape, ComplexMagnetization::zero);

        // Expand tau along mu
        this->_tau.reshape(origin, shape, NAN);

        // Expand p along mu
        Index p_origin(this->_p.origin());
        std::copy(origin.begin(), origin.end(), p_origin.begin()+1);

        Shape p_shape(this->_p.shape());
        std::copy(shape.begin(), shape.end(), p_shape.begin()+1);
        p_shape[0] = 3;

        this->_p.reshape(p_origin, p_shape, NAN);

        // Expand F along mu
        Index F_origin(this->_F.origin());
        std::copy(origin.begin(), origin.end(), F_origin.begin()+1);

        Shape F_shape(this->_F.shape());
        std::copy(shape.begin(), shape.end(), F_shape.begin()+1);
        F_shape[0] = 3;
        F_shape[F_shape.size()-1] = this->_dimensions.size();

        this->_F.reshape(F_origin, F_shape, NAN);
    }

    // Relaxation effects
    auto const E_1 = std::exp(
        (-this->_species.get_R1() * time_interval.get_duration()).magnitude);
    auto const E_2 = std::exp(
        (-this->_species.get_R2() * time_interval.get_duration()).magnitude);

    // Configuration-independent diffusion effects
    Array<Real> p_mu(time_interval.get_gradient_moment().size());
    for(std::size_t i=0; i<p_mu.size(); ++i)
    {
        p_mu[i] = time_interval.get_gradient_moment()[i].convert_to(rad/m);
    }
    auto const p_mu_norm_third = dot(p_mu, p_mu)/3;
    auto const minus_D_tau =
        (-this->_species.get_D()[0]*time_interval.get_duration()).convert_to(m*m);
    auto const has_gradient = p_mu != Array<Real>(p_mu.size(), 0.);
    auto const do_diffusion = this->_species.get_D()[0] > 0*m*m/s && has_gradient;

    // Offsets on the mu axis
    auto && m_stride = this->_m.stride()[mu];
    auto && tau_stride = this->_tau.stride()[mu];
    auto && p_stride = this->_p.stride()[1+mu];
    auto && F_stride = this->_F.stride()[1+mu];

    // Each line along the mu axis can be updated independently. We scan the
    // first hyperplane of the bounding box orthogonal to the mu axis, and for
    // each point of this hyperplane, we scan the mu-line forward (to compute
    // m.m and m.z on each point) and backward (to compute m.p on each point)
    Shape hyperplane_shape(this->_bounding_box.second);
    hyperplane_shape[mu] = 1;
    GridScanner const scanner(
        this->_m.origin(), this->_m.shape(),
        this->_bounding_box.first, hyperplane_shape);
    // NOTE: since GridScanner is not a RandomAccessIterator, we need to
    // pre-compute the values in order to use OpenMP.
    std::vector<GridScanner::value_type> const scanner_data(
        scanner.begin(), scanner.end());
    #pragma omp parallel for
#ifdef _WIN32
    // WARNING: only signed integer types in OpenMP loops on Windows
    for(int scanner_index=0; scanner_index<scanner_data.size(); ++scanner_index)
    {
        auto scanner_it = scanner_data.data()+scanner_index;
#else
    for(
        auto scanner_it=scanner_data.begin(); scanner_it<scanner_data.end();
        ++scanner_it)
    {
#endif
        // Position of the first point of the line
        auto && line_start_index = scanner_it->first;
        auto && line_start_offset = scanner_it->second;

        // Iterator pointing to the first point of the line
        auto m_line_start_it = this->_m.data() + line_start_offset;
        auto tau_line_start_it = this->_tau.data() + line_start_offset;
        auto p_line_start_it = this->_p.data() + 3*line_start_offset;
        auto F_line_start_it =
            this->_F.data()
            + mu*this->_F.stride()[this->_F.dimension()-1]
            + 3*line_start_offset;

        // Update tau and p along the line. F must to be updated during the
        // forward/backward pass since it is updated based on neighbors.
        {
            auto index=line_start_index;
            auto tau_it = tau_line_start_it;
            auto p_it = p_line_start_it;
            for(int i=0; i<this->_bounding_box.second[mu]; ++i)
            {
                if(std::isnan(*tau_it))
                {
                    this->_compute_tau_n(index, *tau_it);
                }
                if(has_gradient)
                {
                    Array<Real> p_n(p_it, 3);
                    if(std::isnan(p_n[0]))
                    {
                        this->_compute_p_n(index, p_n);
                    }
                }
                ++index[mu];
                tau_it += tau_stride;
                p_it += p_stride;
            }
        }

        // Forward scan: compute m.m and m.z
        for(Index::value_type i=0; i<this->_bounding_box.second[mu]-1; ++i)
        {
            // Iterator pointing to the current point of the line
            auto const m_it = m_line_start_it + i*m_stride;
            auto const F_it = F_line_start_it + i*F_stride;

            Real F_minus = 1;
            Real F_zero = 1;

            if(do_diffusion)
            {
                {
                    Array<Real> p_n(p_line_start_it + (i+1)*p_stride, 3);

                    auto const F_minus_it(F_it+F_stride);
                    if(std::isnan(*F_minus_it))
                    {
                        *F_minus_it = this->_compute_F(
                            p_n, -1, p_mu, minus_D_tau, p_mu_norm_third);
                    }
                    F_minus = *F_minus_it;
                }
                {
                    Array<Real> p_n(p_line_start_it + i*p_stride, 3);

                    auto const F_zero_it(F_it+1);
                    if(std::isnan(*F_zero_it))
                    {
                        *F_zero_it = this->_compute_F(
                            p_n, 0, p_mu, minus_D_tau, p_mu_norm_third);
                    }
                    F_zero = *F_zero_it;
                }
            }

            m_it->m = F_minus * E_2 * (m_it+m_stride)->m;
            m_it->z = F_zero * E_1 * m_it->z;
        }

        // Backward scan: compute m.p
        for(Index::value_type i=this->_bounding_box.second[mu]-1; i>0; --i)
        {
            // Iterator pointing to the current point of the line
            auto const m_it = m_line_start_it + i*m_stride;
            auto const F_it = F_line_start_it + i*F_stride;

            Real F_plus = 1;
            if(do_diffusion)
            {
                Array<Real> p_n(p_line_start_it + (i-1)*p_stride, 3);
                auto const F_plus_it(F_it-F_stride+2);
                if(std::isnan(*F_plus_it))
                {
                    *F_plus_it = this->_compute_F(
                        p_n, +1, p_mu, minus_D_tau, p_mu_norm_third);
                }
                F_plus = *F_plus_it;
            }
            m_it->p = F_plus * E_2 * (m_it-m_stride)->p;
        }
    }

    // Repolarization: second term of Equation 19
    auto const repolarization = this->_species.w * (1-E_1);
    this->_m[Index(this->_time_intervals.size(), 0)].p +=
        repolarization * this->_initial_magnetization.p;
    this->_m[Index(this->_time_intervals.size(), 0)].z +=
        repolarization * this->_initial_magnetization.z;
    this->_m[Index(this->_time_intervals.size(), 0)].m +=
        repolarization * this->_initial_magnetization.m;

    this->_timers["time_interval"] +=
        std::chrono::duration<double, std::ratio<1>>(
            std::chrono::high_resolution_clock::now()-start).count();

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
::isochromat(
    std::set<Index> const & configurations, Point const & position,
    Quantity const & relative_frequency) const
{
    auto const start = std::chrono::high_resolution_clock::now();

    Array<Complex> isochromat{0, 0, 0};

    if(position.size() != 3 && !position.empty())
    {
        throw std::runtime_error("Invalid position");
    }
    Array<Real> position_real(position.size(), 0);
    std::transform(
        position.begin(), position.end(), position_real.begin(),
        [](Point::value_type const & x){ return x.convert_to(units::m); });

    auto const omega =
        relative_frequency.convert_to(units::rad/units::s)
        +this->_species.get_delta_omega().convert_to(units::rad/units::s);

    auto const update_isochromat = [&](std::size_t const & offset) {
        auto && tau = this->_tau[offset];

        Real const susceptibility =
            -this->_species.get_R2_prime().convert_to(units::Hz) * tau;

        Real const off_resonance = omega * tau;

        Real gradients_dephasing = 0;
        if(!position_real.empty() && !std::isnan(this->_p[3*offset]))
        {
            Array<Real> const p(const_cast<Real*>(this->_p.data())+3*offset, 3);
            gradients_dephasing = dot(p, position_real);
        }

        auto const dephasing =
            std::exp(Complex(susceptibility, off_resonance-gradients_dephasing));

        auto && m = this->_m[offset];
        isochromat[0] += m.p*dephasing;
        isochromat[1] += m.z*dephasing;
        isochromat[2] += m.m*dephasing;
    };
    auto const update_isochromat_no_offset = [&](Index const & n) {
        auto const offset = dot(n-this->_m.origin(), this->_m.stride());
        return update_isochromat(offset);
    };

    if(configurations.empty())
    {
        GridScanner const scanner(
            this->_m.origin(), this->_m.shape(),
            this->_bounding_box.first, this->_bounding_box.second);
        for(auto && index_offset: scanner)
        {
            update_isochromat(index_offset.second);
        }
    }
    else
    {
        for(auto && configuration: configurations)
        {
            update_isochromat_no_offset(configuration);
        }
    }

    ComplexMagnetization const m_c(
        isochromat[0], isochromat[1].real(), isochromat[2]);
    auto const m_r = as_real_magnetization(m_c);

    this->_timers["isochromat"] +=
        std::chrono::duration<double, std::ratio<1>>(
            std::chrono::high_resolution_clock::now()-start).count();

    return m_r;
}

std::map<std::string, double> const &
Model
::timers() const
{
    return this->_timers;
}

void
Model
::_compute_tau_n(Index const & n, Real & tau)
{
    tau = 0;
    for(std::size_t d=0; d<this->_dimensions.size(); ++d)
    {
        auto && tau_eta = this->_time_intervals[d].get_duration().convert_to(units::s);
        tau += n[d] * tau_eta;
    }
}

void
Model
::_compute_p_n(Index const & n, Array<Real> & p)
{
    std::fill(p.begin(), p.end(), 0);
    for(std::size_t d=0; d<this->_dimensions.size(); ++d)
    {
        Array<Real> p_eta(this->_time_intervals[d].get_gradient_moment().size());
        for(std::size_t i=0; i<p_eta.size(); ++i)
        {
            p_eta[i] =
                this->_time_intervals[d].get_gradient_moment()[i].convert_to(
                    units::rad/units::m);
        }
        p += n[d] * p_eta;
    }
}

Real
Model
::_compute_F(
    Array<Real> const & p_n, int j, Array<Real> const & p_mu,
    Real minus_D_tau, Real p_mu_norm_third)
{
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
        for(std::size_t d=0; d<this->_m.dimension(); ++d)
        {
            first[d] = std::min(first[d], index[d]);
            last[d] = std::max(last[d], index[d]);
        }
    };

    GridScanner const scanner(
        this->_m.origin(), this->_m.shape(),
        this->_bounding_box.first, this->_bounding_box.second);
    for(auto && index_offset: scanner)
    {
        auto && index = index_offset.first;
        auto && offset = index_offset.second;

        if(index == zero_i)
        {
            update_boundary(index);
            continue;
        }

        auto && m = this->_m[offset];
        auto const magnitude =
            m.p*std::conj(m.p) + m.z*m.z + m.m*std::conj(m.m);
        if(magnitude.real() >= this->_epsilon_squared)
        {
            update_boundary(index);
            continue;
        }

        bool would_create_concavity = false;
        for(std::size_t i=0; i<this->_m.dimension(); ++i)
        {
            int neighbors_count = 0;

            auto neighbor_index = index;
            auto neighbor_offset = offset;

            neighbor_index[i] += 1;
            neighbor_offset += this->_m.stride()[i];
            if(
                neighbor_index[i] < this->_m.origin()[i]+this->_m.shape()[i]
                && this->_m[neighbor_offset] != ComplexMagnetization::zero)
            {
                ++neighbors_count;
            }

            neighbor_index[i] -= 2; // NOTE: go two steps back to skip the current index
            neighbor_offset -= 2*this->_m.stride()[i];
            if(
                neighbor_index[i] >= this->_m.origin()[i]
                && this->_m[neighbor_offset] != ComplexMagnetization::zero)
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

        // Otherwise we can remove the configuration.
        m = ComplexMagnetization::zero;
    }

    // Keep a symmetric bounding box as assumed by apply_time_interval
    Index origin(this->_m.dimension());
    Shape shape(this->_m.dimension());
    for(std::size_t i=0; i<shape.size(); ++i)
    {
        origin[i] = -std::max(std::abs(first[i]), std::abs(last[i]));
        shape[i] = 1+2*std::max(std::abs(first[i]), std::abs(last[i]));
    }
    this->_bounding_box = {origin, shape};
}

}

}
