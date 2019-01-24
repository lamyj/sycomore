#include "Model.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include "sycomore/Grid.h"
#include "sycomore/GridScanner.h"
#include "sycomore/magnetization.h"
#include "sycomore/Pulse.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"

namespace sycomore
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

    // Grid of gradient moments
    Index p_origin(origin.size()+1, 0);
    std::copy(origin.begin(), origin.end(), p_origin.begin()+1);

    Shape p_shape(shape.size()+1, 0);
    std::copy(shape.begin(), shape.end(), p_shape.begin()+1);
    p_shape[0] = 3;

    this->_p = Grid<Real>(p_origin, p_shape, NAN);

    // Grid of diffusion damping factors
    Index F_origin(origin.size()+2, 0);
    std::copy(origin.begin(), origin.end(), F_origin.begin()+1);

    Shape F_shape(shape.size()+2, 0);
    std::copy(shape.begin(), shape.end(), F_shape.begin()+1);
    F_shape[0] = 3;
    F_shape[F_shape.size()-1] = this->_dimensions.size();

    this->_F = Grid<Real>(F_origin, F_shape, NAN);
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
    auto const start = std::chrono::high_resolution_clock::now();

    auto const R_m = pulse.rotation_matrix();
    auto && R = R_m.data();

    #pragma omp parallel for
    for(
        auto it=this->_m.data();
        it<this->_m.data()+this->_m.stride()[this->_m.dimension()]; ++it)
    {
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
::apply_pulses(
    std::vector<Pulse> const & pulses, std::string const & interval_name)
{
    if(pulses.empty())
    {
        return;
    }

    this->apply_pulse(pulses[0]);
    for(size_t i=1; i!=pulses.size(); ++i)
    {
        this->apply_time_interval(interval_name);
        this->apply_pulse(pulses[i]);
    }
}

void
Model
::apply_time_interval(std::string const & name)
{
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
    auto const E_1 = std::exp(-this->_species.R1 * time_interval.duration);
    auto const E_2 = std::exp(-this->_species.R2 * time_interval.duration);

    // Configuration-independent diffusion effects
    auto && p_mu = time_interval.gradient_moment;
    auto const p_mu_norm_third = dot(p_mu, p_mu)/3;
    auto const minus_D_tau = -this->_species.D*time_interval.duration;
    auto const has_gradient = p_mu != Array<Real>(p_mu.size(), 0);
    auto const do_diffusion = this->_species.D > 0 && has_gradient;

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
    for(
        auto scanner_it=scanner_data.begin(); scanner_it<scanner_data.end();
        ++scanner_it)
    {
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

        // Forward scan: compute m.m and m.z
        for(Index::value_type i=0; i<this->_bounding_box.second[mu]-1; ++i)
        {
            // Iterator pointing to the current point of the line
            auto const m_it = m_line_start_it + i*m_stride;
            auto const F_it = F_line_start_it + i*F_stride;

            Real F_minus = 1;
            Real F_zero = 1;
            if(has_gradient)
            {
                {
                    Array<Real> p_n(p_line_start_it + (i+1)*p_stride, 3);
                    if(std::isnan(p_n[0]))
                    {
                        auto index(line_start_index);
                        index[mu] += i+1;

                        this->_compute_tau_n(
                            index, *(tau_line_start_it+(i+1)*tau_stride));
                        this->_compute_p_n(index, p_n);
                    }
                    if(do_diffusion)
                    {
                        auto const F_minus_it(F_it+F_stride);
                        if(std::isnan(*F_minus_it))
                        {
                            *F_minus_it = this->_compute_F(
                                p_n, -1, p_mu, minus_D_tau, p_mu_norm_third);
                        }
                        F_minus = *F_minus_it;
                    }
                }
                {
                    Array<Real> p_n(p_line_start_it + i*p_stride, 3);
                    if(std::isnan(p_n[0]))
                    {
                        auto index(line_start_index);
                        index[mu] += i;

                        this->_compute_tau_n(
                            index, *(tau_line_start_it+i*tau_stride));
                        this->_compute_p_n(index, p_n);
                    }
                    if(do_diffusion)
                    {
                        auto const F_zero_it(F_it+1);
                        if(std::isnan(*F_zero_it))
                        {
                            *F_zero_it = this->_compute_F(
                                p_n, 0, p_mu, minus_D_tau, p_mu_norm_third);
                        }
                        F_zero = *F_zero_it;
                    }
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
            if(has_gradient)
            {
                Array<Real> p_n(p_line_start_it + (i-1)*p_stride, 3);
                if(std::isnan(p_n[0]))
                {
                    auto index(line_start_index);
                    index[mu] += i-1;

                    this->_compute_tau_n(
                        index, *(tau_line_start_it+(i-1)*tau_stride));
                    this->_compute_p_n(index, p_n);
                }
                if(do_diffusion)
                {
                    auto const F_plus_it(F_it-F_stride+2);
                    if(std::isnan(*F_plus_it))
                    {
                        *F_plus_it = this->_compute_F(
                            p_n, +1, p_mu, minus_D_tau, p_mu_norm_third);
                    }
                    F_plus = *F_plus_it;
                }
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
    std::set<Index> const & configurations, Point const & position) const
{
    auto const start = std::chrono::high_resolution_clock::now();

    Magnetization isochromat{0, 0, 0};

    Array<Real> position_real(position.size(), 0);
    std::transform(
        position.begin(), position.end(), position_real.begin(),
        [](Point::value_type const & x){ return x.value; });

    auto const update_isochromat = [&](Index const & n, size_t const & offset) {
        Real off_resonance = 0;
        Real gradients_dephasing= 0;
        if(!position_real.empty())
        {
            Array<Real> const p(const_cast<Real*>(this->_p.data())+3*offset, 3);
            gradients_dephasing = dot(p, position_real);
        }

        auto && m = this->_m[n];
        Magnetization m_r;
        if(off_resonance != 0 || gradients_dephasing != 0)
        {
            auto const dephasing =
                std::exp(Complex(0,1)*(off_resonance-gradients_dephasing));

            // NOTE: the additional dephasing (exponential term in eq. 4) will
            // create complex longitudinal magnetization. However, the summation
            // of the longitudinal magnetization will be purely real.
            ComplexMagnetization const dephased_m(
                m.p*dephasing, (m.z*dephasing).real(), m.m*dephasing);
            m_r = as_real_magnetization(dephased_m);
        }
        else
        {
            m_r = as_real_magnetization(m);
        }

        isochromat.x += m_r.x;
        isochromat.y += m_r.y;
        isochromat.z += m_r.z;
    };
    auto const update_isochromat_no_offset = [&](Index const & n) {
        auto const offset = dot(n-this->_m.origin(), this->_m.stride());
        return update_isochromat(n, offset);
    };

    if(configurations.empty())
    {
        GridScanner const scanner(
            this->_m.origin(), this->_m.shape(),
            this->_bounding_box.first, this->_bounding_box.second);
        for(auto && index_offset: scanner)
        {
            update_isochromat(index_offset.first, index_offset.second);
        }
    }
    else
    {
        for(auto && configuration: configurations)
        {
            update_isochromat_no_offset(configuration);
        }
    }

    this->_timers["isochromat"] +=
        std::chrono::duration<double, std::ratio<1>>(
            std::chrono::high_resolution_clock::now()-start).count();

    return isochromat;
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
    for(size_t d=0; d<this->_dimensions.size(); ++d)
    {
        auto && tau_eta = this->_time_intervals[d].duration;
        tau += tau_eta;
    }
}

void
Model
::_compute_p_n(Index const & n, Array<Real> & p)
{
    std::fill(p.begin(), p.end(), 0);
    for(size_t d=0; d<this->_dimensions.size(); ++d)
    {
        auto && p_eta = this->_time_intervals[d].gradient_moment;
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
        for(size_t d=0; d<this->_m.dimension(); ++d)
        {
            first[d] = std::min(first[d], index[d]);
            last[d] = std::max(first[d], index[d]);
        }
    };

    for(auto && index: GridScanner(this->_bounding_box.first, this->_bounding_box.second))
    {
        if(index.first == zero_i)
        {
            update_boundary(index.first);
            continue;
        }

        auto && m = this->_m[index.first];
        auto const magnitude =
            m.p * std::conj(m.p)
            + m.z*m.z
            + m.m * std::conj(m.m);
        if(magnitude.real() >= this->_epsilon_squared)
        {
            update_boundary(index.first);
            continue;
        }

        bool would_create_concavity = false;
        for(size_t i=0; i<this->_m.dimension(); ++i)
        {
            int neighbors_count = 0;

            Index neighbor = index.first;

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
            update_boundary(index.first);
            continue;
        }

        this->_m[index.first] = ComplexMagnetization::zero;
    }

    // FIXME: make sure we keep a symmetric bounding box.
    Shape shape(first.size());
    std::transform(
        first.begin(), first.end(), last.begin(), shape.begin(),
        [](int f, int l) { return l-f+1;});
    this->_bounding_box = {first, shape};
}

}
