#ifndef _faa6a046_30f6_4e87_91a6_033e2330b405
#define _faa6a046_30f6_4e87_91a6_033e2330b405

#include <tuple>
#include <utility>
#include <vector>

#include "sycomore/Quantity.h"
#include "sycomore/Species.h"

namespace sycomore
{

namespace epg
{

namespace operators
{

std::vector<Complex> pulse(Quantity angle, Quantity phase);

std::pair<Real, Real> relaxation(
    Species const & species, Quantity const & duration);

std::tuple<Real, Real, Real> diffusion(
    Species const & species, Quantity const & duration, Quantity const & k, 
    Quantity const & delta_k);

}
    
}

}

#endif // _faa6a046_30f6_4e87_91a6_033e2330b405
