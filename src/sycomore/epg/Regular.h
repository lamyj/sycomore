#ifndef _fbf381fe_fd75_427e_88de_a033418c943c
#define _fbf381fe_fd75_427e_88de_a033418c943c

#include <vector>

#include "sycomore/magnetization.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace sycomore
{

namespace epg
{

class Regular
{
public:
    Species species;
    Quantity gamma;
    
    Regular(
        Species const & species, 
        Magnetization const & initial_magnetization={0,0,1}, 
        unsigned int initial_size=100, 
        Quantity gamma=2*M_PI*units::rad * 42.57747892*units::MHz/units::T);
    
    Regular(Regular const &) = default;
    Regular(Regular &&) = default;
    Regular & operator=(Regular const &) = default;
    Regular & operator=(Regular &&) = default;
    ~Regular() = default;
    
    std::size_t const states_count() const;
    std::vector<Complex> magnetization(std::size_t state) const;
    Complex const & echo() const;
    
    void apply_pulse(Quantity angle, Quantity phase=0*units::rad);
    void apply_time_interval(
        Quantity const & duration, 
        Quantity const & gradient=0*units::T/units::m);
    void apply_gradient();
    void apply_relaxation(Quantity const & duration);
    void apply_diffusion(Quantity const & duration, Quantity const & gradient);
private:
    std::vector<Complex> _magnetization;
    unsigned int _states_count;
};

}

}

#endif // _fbf381fe_fd75_427e_88de_a033418c943c
