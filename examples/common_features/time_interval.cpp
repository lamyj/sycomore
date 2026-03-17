#include <iostream>

#include <sycomore/sycomore.h>
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

#if __has_include(<xtensor/xtensor.hpp>)
#include <xtensor/xio.hpp>
#else
#include <xtensor/io/xio.hpp>
#endif

int main()
{
    using namespace sycomore::units;
    
    // Scalar gradient, defined by its amplitude
    sycomore::TimeInterval interval(1*ms, 20*mT/m);
    std::cout
        << interval.duration() << "\n"
        << interval.gradient_amplitude() << "\n"
        << interval.gradient_area()/interval.duration() << "\n"
        << interval.gradient_dephasing()/(sycomore::gamma*interval.duration())
        << "\n";
    return 0;
}
