#ifndef _d9169a5f_d53b_4440_bfc7_2b3f978b665d
#define _d9169a5f_d53b_4440_bfc7_2b3f978b665d

#include <vector>

#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

class Discrete
{
public:
    Species species;
    Quantity gamma;
    
    Discrete(
        Species const & species, 
        Magnetization const & initial_magnetization={0,0,1}, 
        Quantity bin_width=1*units::rad/units::m, 
        Quantity gamma=2*M_PI*units::rad * 42.57747892*units::MHz/units::T);
    
    Discrete(Discrete const &) = default;
    Discrete(Discrete &&) = default;
    Discrete & operator=(Discrete const &) = default;
    Discrete & operator=(Discrete &&) = default;
    ~Discrete() = default;
    
    std::vector<int64_t> const & orders() const;
    std::vector<Complex> magnetization(std::size_t state_index) const;
    Complex const & echo() const;
    
    void apply_pulse(Quantity angle, Quantity phase=0*units::rad);
    void apply_time_interval(
        Quantity const & duration, 
        Quantity const & gradient=0*units::T/units::m, Real threshold=0);
    void apply_gradient(Quantity const & duration, Quantity const & gradient);
    void apply_relaxation(Quantity const & duration);
    void apply_diffusion(Quantity const & duration, Quantity const & gradient);

private:
    std::vector<Complex> _magnetization;
    Quantity _bin_width;
    std::vector<int64_t> _orders;
};
    
}

}

#endif // _d9169a5f_d53b_4440_bfc7_2b3f978b665d
