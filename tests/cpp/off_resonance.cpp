#define BOOST_TEST_MODULE OffResonance
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>

#include <sycomore/magnetization.h>
#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

using namespace sycomore::units;

BOOST_AUTO_TEST_CASE(OffResonance, *boost::unit_test::tolerance(1e-12))
{
    sycomore::Species const species(0_Hz, 0_Hz, 0_um*um/ms);
    sycomore::Magnetization const m0{0,0,1};

    sycomore::Pulse const pulse(90_deg, M_PI*rad);
    auto const pulse_duration=1_ms;
    int const pulse_support_size = 100;
    int const zero_crossings = 2;

    // NOTE: in the absence of relaxation and diffusion, the TR is meaningless
    auto const TR=500_ms;
    auto const slice_thickness=1_mm;

    auto const pulses = sycomore::hard_pulse_approximation(
        pulse,
        [](sycomore::Real x) { return x==0?1:std::sin(x*M_PI)/(x*M_PI); },
        sycomore::linspace(
            -sycomore::Real(zero_crossings), +sycomore::Real(zero_crossings),
            1+pulse_support_size));

    auto const bandwidth = 2 * zero_crossings / pulse_duration;
    auto const slice_selection_gradient_moment =
        2*M_PI*bandwidth*pulse_duration/slice_thickness;

    auto const frequencies = sycomore::linspace(-30._rad/ms, 30._rad/ms, 201);

    sycomore::Model model(
        species, m0, {
            {"rf", {
                pulse_duration/pulse_support_size, {
                    0*1/m, 0*1/m,
                    slice_selection_gradient_moment/pulse_support_size}}},
            {"refocalization", {
                (TR-pulse_duration)/2., {
                    0*1/m, 0*1/m,
                    -slice_selection_gradient_moment/2}}}
    });

    model.apply_pulses(pulses, "rf");
    model.apply_time_interval("refocalization");

    std::vector<sycomore::Magnetization> magnetization;
    for(auto && frequency: frequencies)
    {
        auto const signal = model.isochromat({}, {}, frequency);
        magnetization.push_back(signal);
    }

    std::vector<double> baseline;
    std::string const root(getenv("SYCOMORE_TEST_DATA")?getenv("SYCOMORE_TEST_DATA"):"");
    if(root.empty())
    {
        throw std::runtime_error("SYCOMORE_TEST_DATA is undefined");
    }
    std::ifstream stream(root+"/baseline/off_resonance.dat", std::ios_base::binary);
    while(stream.good())
    {
        double value;
        stream.read(reinterpret_cast<char*>(&value), sizeof(value));
        if(stream.good())
        {
            baseline.push_back(value);
        }
    }

    // WARNING: we are using absolute tolerance, not relative to the value of
    // left and right
    BOOST_REQUIRE_EQUAL(baseline.size(), 2*magnetization.size());
    for(size_t i=0; i<magnetization.size(); ++i)
    {
        BOOST_TEST(
            magnetization[i].transversal()-baseline[2*i] == 0.,
            "Error at " << i << " [ "
                << magnetization[i].transversal() << " != " << baseline[2*i] << " ]");
        BOOST_TEST(
            magnetization[i].z-baseline[2*i+1] == 0.,
            "Error at " << i << " [ "
                << magnetization[i].z << " != " << baseline[2*i+1] << " ]");
    }
}
