#define BOOST_TEST_MODULE GRE
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <fstream>

#include <sycomore/magnetization.h>
#include <sycomore/como/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

#include <iostream>

using namespace sycomore::units;

struct Fixture
{
    static sycomore::Species const species;
    static sycomore::Magnetization const m0;

    static sycomore::Quantity const flip_angle;
    static sycomore::Quantity const pulse_duration;
    static int const pulse_support_size = 101;
    static int const zero_crossings = 2;

    static sycomore::Quantity const TR;
    static sycomore::Quantity const slice_thickness;

    static size_t const TR_count = 10;
};
sycomore::Species const Fixture::species{1000_ms, 100_ms, 0.89_um*um/ms};
sycomore::Magnetization const Fixture::m0{0,0,1};
sycomore::Quantity const Fixture::flip_angle=40_deg;
sycomore::Quantity const Fixture::pulse_duration=1_ms;
sycomore::Quantity const Fixture::TR=500_ms;
sycomore::Quantity const Fixture::slice_thickness=1_mm;

BOOST_FIXTURE_TEST_CASE(Ideal, Fixture, *boost::unit_test::tolerance(1e-9))
{
    sycomore::como::Model model(species, m0, {{"echo", {TR/2.}}});

    std::vector<sycomore::Magnetization> magnetization;
    for(int i=0; i<TR_count; ++i)
    {
        model.apply_pulse({flip_angle, (M_PI/3+(i%2)*M_PI)*rad});
        model.apply_time_interval("echo");
        magnetization.push_back(model.isochromat());
        model.apply_time_interval("echo");
    }

    std::vector<double> baseline;
    std::string const root(getenv("SYCOMORE_TEST_DATA")?getenv("SYCOMORE_TEST_DATA"):"");
    if(root.empty())
    {
        throw std::runtime_error("SYCOMORE_TEST_DATA is undefined");
    }
    std::ifstream stream(root+"/baseline/GRE_ideal.dat", std::ios_base::binary);
    while(stream.good())
    {
        double value;
        stream.read(reinterpret_cast<char*>(&value), sizeof(value));
        if(stream.good())
        {
            baseline.push_back(value);
        }
    }

    BOOST_REQUIRE_EQUAL(baseline.size(), 3*TR_count);
    for(size_t i=0; i!=TR_count; ++i)
    {
        BOOST_TEST(magnetization[i][0] == baseline[3*i+0]);
        BOOST_TEST(magnetization[i][1] == baseline[3*i+1]);
        BOOST_TEST(magnetization[i][2] == baseline[3*i+2]);
    }
}

BOOST_FIXTURE_TEST_CASE(Real, Fixture, *boost::unit_test::tolerance(1e-9))
{
    auto const t0 = pulse_duration/(2*zero_crossings);
    sycomore::HardPulseApproximation sinc_pulse(
        {flip_angle, 0_rad},
        sycomore::linspace(pulse_duration, pulse_support_size),
        sycomore::sinc_envelope(t0), 1/t0, slice_thickness, "rf");

    sycomore::como::Model model(
        species, m0, {
            {"rf", sinc_pulse.get_time_interval()},
            {"half_echo", {
                (TR-pulse_duration)/2., -sinc_pulse.get_gradient_moment()/2}}
    });

    std::vector<sycomore::Magnetization> magnetization;

    // CoMoTk: 6.13 s
    // Total time: 0.85933 s
    // Pulse time: 0.334359 s
    // Time interval time: 0.521324 s
    // Unaccounted time: 0.00364655 s
    // Sycomore speedup: 7.13


    auto const start = std::chrono::high_resolution_clock::now();
    for(int i=0; i<TR_count; ++i)
    {
        sinc_pulse.set_phase((M_PI/3+(i%2)*M_PI)*rad);
        model.apply_pulse(sinc_pulse);
        model.apply_time_interval("half_echo");
        magnetization.push_back(model.isochromat());
        model.apply_time_interval("half_echo");
    }
    auto const total_time = std::chrono::duration<double, std::ratio<1>>(
        std::chrono::high_resolution_clock::now()-start).count();

    std::cout << "{\"total\": " << total_time;
    double accounted=0;
    for(auto && item: model.timers())
    {
        std::cout << ", \"" << item.first << "\": " << item.second;
        accounted += item.second;
    }
    std::cout << ", \"extra\": " << total_time-accounted << "}\n";

    std::vector<double> baseline;
    std::string const root(getenv("SYCOMORE_TEST_DATA")?getenv("SYCOMORE_TEST_DATA"):"");
    if(root.empty())
    {
        throw std::runtime_error("SYCOMORE_TEST_DATA is undefined");
    }
    std::ifstream stream(root+"/baseline/GRE_real.dat", std::ios_base::binary);
    while(stream.good())
    {
        double value;
        stream.read(reinterpret_cast<char*>(&value), sizeof(value));
        if(stream.good())
        {
            baseline.push_back(value);
        }
    }

    BOOST_REQUIRE_EQUAL(baseline.size(), 3*TR_count);
    for(size_t i=0; i!=TR_count; ++i)
    {
        BOOST_TEST(magnetization[i][0] == baseline[3*i+0]);
        BOOST_TEST(magnetization[i][1] == baseline[3*i+1]);
        BOOST_TEST(magnetization[i][2] == baseline[3*i+2]);
    }
}
