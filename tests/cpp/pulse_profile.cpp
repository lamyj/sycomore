#define BOOST_TEST_MODULE PulseProfile
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <fstream>

#include <sycomore/magnetization.h>
#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

using namespace sycomore::units;

BOOST_AUTO_TEST_CASE(PulseProfile, *boost::unit_test::tolerance(1e-10))
{
    // FIXME: need value of D to compute p_n
    sycomore::Species const species(1000_ms, 100_ms, 0.89_um*um/ms);
    sycomore::Magnetization const m0{0,0,1};

    sycomore::units::Angle const flip_angle=90_deg;
    sycomore::units::Time const pulse_duration=1_ms;
    int const pulse_support_size = 100;
    int const zero_crossings = 2;

    sycomore::units::Time const TR=500_ms;
    sycomore::units::Length const slice_thickness=1_mm;

    int const sampling_support_size = 501;

    sycomore::Array<sycomore::Real> small_flip_angles(pulse_support_size+1, 0.);
    for(size_t i=0; i<small_flip_angles.size(); ++i)
    {
        auto const x = zero_crossings * (-1. + 2.*i/pulse_support_size);
        auto const y = (x==0?1:std::sin(x*M_PI)/(x*M_PI));
        small_flip_angles[i] = y;
    }
    small_flip_angles *=
        flip_angle.convert_to(rad)
        / std::accumulate(small_flip_angles.begin(), small_flip_angles.end(), 0.);

    auto const bandwidth = 2 * zero_crossings / pulse_duration;
    auto const slice_selection_gradient_moment =
        2*M_PI*bandwidth*pulse_duration/slice_thickness;

    std::vector<sycomore::Array<sycomore::Real>> sampling_locations;
    for(size_t i=0; i<sampling_support_size; ++i)
    {
        auto const delta = 2.*slice_thickness/(sampling_support_size-1);
        auto const z = -slice_thickness + i*delta;
        sampling_locations.push_back({0, 0, z.convert_to(m)});
    }

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

    auto const phase = M_PI*rad;

    model.apply_pulse({small_flip_angles[0]*rad, phase});
    for(size_t j=1; j!=small_flip_angles.size(); ++j)
    {
        model.apply_time_interval("rf");
        model.apply_pulse({small_flip_angles[j]*rad, phase});
    }

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

#define TEST_COMPONENT(left, right, where) \
    BOOST_TEST(\
        left == right, \
        "Error on " << #left << " (" << where << ") at " << x \
            << " [ " << left << " != " << right << " ]")
#define TEST_MAGNETIZATION(where) \
    TEST_COMPONENT(m.x, *(baseline_it+0), where); \
    TEST_COMPONENT(m.y, *(baseline_it+1), where); \
    TEST_COMPONENT(m.z, *(baseline_it+2), where)

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
