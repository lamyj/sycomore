#define BOOST_TEST_MODULE PulseProfile
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>

#include <sycomore/HardPulseApproximation.h>
#include <sycomore/magnetization.h>
#include <sycomore/como/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

using namespace sycomore::units;

BOOST_AUTO_TEST_CASE(PulseProfile, *boost::unit_test::tolerance(1e-9))
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

    int const sampling_support_size = 501;

    auto const t0 = pulse_duration/(2*zero_crossings);
    sycomore::HardPulseApproximation const sinc_pulse(
        pulse,
        sycomore::linspace(pulse_duration, pulse_support_size),
        sycomore::sinc_envelope(t0), 1/t0, slice_thickness, "rf");

    sycomore::TimeInterval const refocalization(
        (TR-pulse_duration)/2., -sinc_pulse.get_gradient_moment()/2);

    auto const sampling_locations = sycomore::linspace(
        sycomore::Point{0_m, 0_m, 2*slice_thickness},
        sampling_support_size);

    sycomore::como::Model model(
        species, m0, {
            {"rf", sinc_pulse.get_time_interval()},
            {"refocalization", refocalization}
    });

    model.apply_pulse(sinc_pulse);

    std::vector<sycomore::Magnetization> before_refocalization;
    for(auto && location: sampling_locations)
    {
        auto const signal = model.isochromat({}, location);
        before_refocalization.push_back(signal);
    }

    model.apply_time_interval("refocalization");

    std::vector<sycomore::Magnetization> after_refocalization;
    for(auto && location: sampling_locations)
    {
        auto const signal = model.isochromat({}, location);
        after_refocalization.push_back(signal);
    }

    std::vector<double> baseline;
    std::string const root(getenv("SYCOMORE_TEST_DATA")?getenv("SYCOMORE_TEST_DATA"):"");
    if(root.empty())
    {
        throw std::runtime_error("SYCOMORE_TEST_DATA is undefined");
    }
    std::ifstream stream(root+"/baseline/pulse_profile.dat", std::ios_base::binary);
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
#define TEST_COMPONENT(left, right, where) \
    BOOST_TEST(\
    left-right == 0., \
    "Error on " << #left << " (" << where << ") at " << x \
        << " [ " << left << " != " << right << " ]")
#define TEST_MAGNETIZATION(where) \
    TEST_COMPONENT(m[0], *(baseline_it+0), where); \
    TEST_COMPONENT(m[1], *(baseline_it+1), where); \
    TEST_COMPONENT(m[2], *(baseline_it+2), where)

    BOOST_REQUIRE_EQUAL(baseline.size(), 2*3*sampling_locations.size());
    auto baseline_it = baseline.begin();
    for(auto && m: before_refocalization)
    {
        auto x = sampling_locations[(baseline_it-baseline.begin())/3];
        TEST_MAGNETIZATION("before");
        baseline_it += 3;
    }
    for(auto && m: after_refocalization)
    {
        auto x = sampling_locations[
            (baseline_it-baseline.begin())/3-before_refocalization.size()];
        TEST_MAGNETIZATION("after");
        baseline_it += 3;
    }
}
