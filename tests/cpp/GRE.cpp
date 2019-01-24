#define BOOST_TEST_MODULE GRE
#include <boost/test/unit_test.hpp>

#include <chrono>
#include <cmath>
#include <fstream>

#include <sycomore/magnetization.h>
#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

#include <iostream>

using namespace sycomore::units;

struct Fixture
{
    static sycomore::Species const species;
    static sycomore::Magnetization const m0;

    static sycomore::units::Angle const flip_angle;
    static sycomore::units::Time const pulse_duration;
    static int const pulse_support_size = 100;
    static int const zero_crossings = 2;

    static sycomore::units::Time const TR;
    static sycomore::units::Length const slice_thickness;

    static size_t const TR_count = 10;
};
sycomore::Species const Fixture::species{1000_ms, 100_ms, 0.89_um*um/ms};
sycomore::Magnetization const Fixture::m0{0,0,1};
sycomore::units::Angle const Fixture::flip_angle=40_deg;
sycomore::units::Time const Fixture::pulse_duration=1_ms;
sycomore::units::Time const Fixture::TR=500_ms;
sycomore::units::Length const Fixture::slice_thickness=1_mm;

BOOST_FIXTURE_TEST_CASE(Ideal, Fixture, *boost::unit_test::tolerance(1e-15))
{
    sycomore::Model model(species, m0, {{"echo", {TR/2.}}});

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
        BOOST_TEST(magnetization[i].x == baseline[3*i+0]);
        BOOST_TEST(magnetization[i].y == baseline[3*i+1]);
        BOOST_TEST(magnetization[i].z == baseline[3*i+2]);
    }
}

BOOST_FIXTURE_TEST_CASE(Real, Fixture, *boost::unit_test::tolerance(1e-14))
{
    auto const sinc = [](sycomore::Real x) { return x==0?1:std::sin(x*M_PI)/(x*M_PI); };
    auto const pulse_support = sycomore::linspace(
        -sycomore::Real(zero_crossings), +sycomore::Real(zero_crossings),
        1+pulse_support_size);

    auto const bandwidth = 2 * zero_crossings / pulse_duration;
    sycomore::GradientMoment const slice_selection_gradient_moment =
        2*M_PI*bandwidth*pulse_duration/slice_thickness;

    sycomore::Model model(
        species, m0, {
            {"rf", {
                pulse_duration/pulse_support_size, {
                    sycomore::GradientMoment(0), sycomore::GradientMoment(0),
                    slice_selection_gradient_moment/pulse_support_size}}},
            {"echo", {
                (TR-pulse_duration)/2., {
                    sycomore::GradientMoment(0), sycomore::GradientMoment(0),
                    -slice_selection_gradient_moment/2}}}
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
        auto const phase = (M_PI/3+(i%2)*M_PI)*rad;
        auto const pulses = sycomore::hard_pulse_approximation(
            {flip_angle, phase}, sinc, pulse_support);
        model.apply_pulses(pulses, "rf");
        model.apply_time_interval("echo");
        magnetization.push_back(model.isochromat());
        model.apply_time_interval("echo");
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
        BOOST_TEST(magnetization[i].x == baseline[3*i+0]);
        BOOST_TEST(magnetization[i].y == baseline[3*i+1]);
        BOOST_TEST(magnetization[i].z == baseline[3*i+2]);
    }
}
