#define BOOST_TEST_MODULE OffResonance
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>

#include <sycomore/magnetization.h>
#include <sycomore/como/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

using namespace sycomore::units;

BOOST_AUTO_TEST_CASE(OffResonance, *boost::unit_test::tolerance(1e-9))
{
    sycomore::Species const species(0_Hz, 0_Hz, 0_um*um/ms);
    sycomore::Magnetization const m0{0,0,1};

    sycomore::Pulse const pulse(90_deg, M_PI*rad);
    auto const pulse_duration=1_ms;
    int const pulse_support_size = 101;
    int const zero_crossings = 2;

    // NOTE: in the absence of relaxation and diffusion, the TR is meaningless
    auto const TR=500_ms;
    auto const slice_thickness=1_mm;

    auto const t0 = pulse_duration/(2*zero_crossings);
    sycomore::HardPulseApproximation const sinc_pulse(
        pulse,
        sycomore::linspace(pulse_duration, pulse_support_size),
        sycomore::sinc_envelope(t0), 1/t0, slice_thickness, "rf");

    auto const frequencies = sycomore::linspace(60._rad/ms, 201);

    sycomore::como::Model model(
        species, m0, {
            {"rf", sinc_pulse.get_time_interval()},
            {"refocalization", {
                (TR-pulse_duration)/2., -sinc_pulse.get_gradient_moment()/2}}
    });

    model.apply_pulse(sinc_pulse);
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
            sycomore::transversal(magnetization[i])-baseline[2*i] == 0.,
            "Error at " << i << " [ "
                << sycomore::transversal(magnetization[i]) <<
                " != " << baseline[2*i] << " ]");
        BOOST_TEST(
            magnetization[i][2]-baseline[2*i+1] == 0.,
            "Error at " << i << " [ "
                << magnetization[i][2] << " != " << baseline[2*i+1] << " ]");
    }
}
