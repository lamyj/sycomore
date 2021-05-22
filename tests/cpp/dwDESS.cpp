#define BOOST_TEST_MODULE dwDESS
#include <boost/test/unit_test.hpp>

#include <utility>

#include <sycomore/epg/Discrete.h>
#include <sycomore/epg/Regular.h>
#include <sycomore/Quantity.h>
#include <sycomore/Species.h>
#include <sycomore/units.h>


#include <boost/math/special_functions/relative_difference.hpp>

using namespace sycomore::units;

/*
 * Return the DESS signal with treatment of diffusion, cf. Freed et al., 
 * Steady-state free precession experiments and exact treatment of diffusion in
 * a uniform gradient. Journal of Chemical Physics 115(9), pp. 4249-4258. 2001.
 * doi:10.1063/1.1389859
 */
std::pair<double, double> freed(
    sycomore::Species const & species, sycomore::Quantity const & alpha, 
    sycomore::Quantity const & G_tau, sycomore::Quantity const & TR,
    double M0=1, int truncation_level=6);

struct Fixture
{
    // This follows the DW-DESS parameters in Materials&Methods/Numerical 
    // simulations of "Quantitative In Vivo Diffusion Imaging of Cartilage
    // Using Double Echo Steady-State Free Precession", Bieri & al., Magnetic
    // Resonance in Medicine 68(3), 2012.
    // doi:10.1002/mrm.23275
    
    sycomore::Species const species{1200_ms, 50_ms, 1_um*um/ms};

    sycomore::Quantity const flip_angle{20_deg};
    sycomore::Quantity const G_tau{200_mT/m*ms};
    sycomore::Quantity const TR{10_ms};
        
    std::pair<double, double> const analytical{
        freed(species, flip_angle, G_tau, TR)};
};

BOOST_FIXTURE_TEST_CASE(Regular, Fixture, *boost::unit_test::tolerance(1e-5))
{
    sycomore::epg::Regular model(species, {0,0,1}, 100, G_tau);
    std::pair<double, double> epg;
    unsigned int const repetitions = 5*species.get_T1()/TR;
    for(unsigned int repetition = 0; repetition != repetitions; ++repetition)
    {
        model.apply_pulse(flip_angle);
        epg.first = std::abs(model.echo());
        model.apply_time_interval(TR, G_tau/TR);
        epg.second = std::abs(model.echo());
    }
    
    BOOST_TEST(analytical.first == epg.first);
    BOOST_TEST(analytical.second == epg.second);
}

BOOST_FIXTURE_TEST_CASE(RegularThreshold, Fixture, *boost::unit_test::tolerance(1e-5))
{
    sycomore::epg::Regular model(species, {0,0,1}, 100, G_tau);
    model.threshold = 1e-6;
    
    std::pair<double, double> epg;
    unsigned int const repetitions = 5*species.get_T1()/TR;
    for(unsigned int repetition = 0; repetition != repetitions; ++repetition)
    {
        model.apply_pulse(flip_angle);
        epg.first = std::abs(model.echo());
        model.apply_time_interval(TR, G_tau/TR);
        epg.second = std::abs(model.echo());
    }
    
    BOOST_TEST(analytical.first == epg.first);
    BOOST_TEST(analytical.second == epg.second);
}

BOOST_FIXTURE_TEST_CASE(DiscreteNoReadout, Fixture, *boost::unit_test::tolerance(1e-5))
{
    sycomore::epg::Discrete model(species);
    std::pair<double, double> epg;
    unsigned int const repetitions = 5*species.get_T1()/TR;
    for(unsigned int repetition = 0; repetition != repetitions; ++repetition)
    {
        model.apply_pulse(flip_angle);
        epg.first = std::abs(model.echo());
        model.apply_time_interval(TR, G_tau/TR);
        epg.second = std::abs(model.echo());
    }
    
    BOOST_TEST(analytical.first == epg.first);
    BOOST_TEST(analytical.second == epg.second);
}

BOOST_FIXTURE_TEST_CASE(DiscreteThresholdNoReadout, Fixture, *boost::unit_test::tolerance(1e-5))
{
    sycomore::epg::Discrete model(species);
    model.threshold = 1e-6;
    
    std::pair<double, double> epg;
    unsigned int const repetitions = 5*species.get_T1()/TR;
    for(unsigned int repetition = 0; repetition != repetitions; ++repetition)
    {
        model.apply_pulse(flip_angle);
        epg.first = std::abs(model.echo());
        model.apply_time_interval(TR, G_tau/TR);
        epg.second = std::abs(model.echo());
    }
    
    BOOST_TEST(analytical.first == epg.first);
    BOOST_TEST(analytical.second == epg.second);
}

BOOST_FIXTURE_TEST_CASE(DiscreteWithReadout, Fixture, *boost::unit_test::tolerance(5e-5))
{
    // This follows the DW-DESS implementation in Fig. 2 of "Quantitative In 
    // Vivo Diffusion Imaging of Cartilage Using Double Echo Steady-State Free
    // Precession", Bieri & al., Magnetic Resonance in Medicine 68(3), 2012.
    // doi:10.1002/mrm.23275
    // The signal will be be affected by
    // - The shift from start and end of TR (readout duration)
    // - The readout gradient area
    // Using unrealistic values for readout duration and voxel size as well as
    // correcting the diffusion time interval to include readout gradient allows
    // staying close to the analytical model.
    
    sycomore::Quantity const readout_duration = 0.1_us;
    sycomore::Quantity const voxel_size{10_cm};
    sycomore::TimeInterval const ro_plus{
        0.5*readout_duration, 0.5 * 2*M_PI/voxel_size};
    sycomore::TimeInterval const ro_minus{
        ro_plus.get_duration(), -ro_plus.get_gradient_amplitude()};
    sycomore::TimeInterval const diffusion{
        TR-6*readout_duration, G_tau - 2*ro_plus.get_gradient_area()[0]};

    sycomore::epg::Discrete model(species);
    std::pair<double, double> epg;
    unsigned int const repetitions = 5*species.get_T1()/TR;
    for(unsigned int repetition = 0; repetition != repetitions; ++repetition)
    {
        model.apply_pulse(flip_angle);

        model.apply_time_interval(ro_minus);
        model.apply_time_interval(ro_plus);
        epg.first = std::abs(model.echo());
        model.apply_time_interval(ro_plus);

        model.apply_time_interval(diffusion);

        model.apply_time_interval(ro_plus);
        epg.second = std::abs(model.echo());
        model.apply_time_interval(ro_plus);
        model.apply_time_interval(ro_minus);
    }

    BOOST_TEST(analytical.first == epg.first);
    BOOST_TEST(analytical.second == epg.second);
}

BOOST_FIXTURE_TEST_CASE(DiscreteThresholdWithReadout, Fixture, *boost::unit_test::tolerance(5e-5))
{
    sycomore::Quantity const readout_duration = 0.1_us;
    sycomore::Quantity const voxel_size{10_cm};
    sycomore::TimeInterval const ro_plus{
        0.5*readout_duration, 0.5 * 2*M_PI/voxel_size};
    sycomore::TimeInterval const ro_minus{
        ro_plus.get_duration(), -ro_plus.get_gradient_amplitude()};
    sycomore::TimeInterval const diffusion{
        TR-6*readout_duration, G_tau - 2*ro_plus.get_gradient_area()[0]};

    sycomore::epg::Discrete model(species);
    model.threshold = 1e-6;
    
    std::pair<double, double> epg;
    unsigned int const repetitions = 5*species.get_T1()/TR;
    for(unsigned int repetition = 0; repetition != repetitions; ++repetition)
    {
        model.apply_pulse(flip_angle);

        model.apply_time_interval(ro_minus);
        model.apply_time_interval(ro_plus);
        epg.first = std::abs(model.echo());
        model.apply_time_interval(ro_plus);

        model.apply_time_interval(diffusion);

        model.apply_time_interval(ro_plus);
        epg.second = std::abs(model.echo());
        model.apply_time_interval(ro_plus);
        model.apply_time_interval(ro_minus);
    }

    BOOST_TEST(analytical.first == epg.first);
    BOOST_TEST(analytical.second == epg.second);
}

std::pair<double, double> freed(
    sycomore::Species const & species, sycomore::Quantity const & alpha, 
    sycomore::Quantity const & G_tau, sycomore::Quantity const & TR_,
    double M0, int truncation_level)
{
    using namespace sycomore::units;
    
    // Cast early to avoid repetitions in the continued fraction
    auto const R1 = species.get_R1().convert_to(Hz);
    auto const R2 = species.get_R2().convert_to(Hz);
    auto const TR = TR_.convert_to(s);
    
    auto const cos_alpha = std::cos(alpha.convert_to(rad));
    auto const sin_alpha = std::sin(alpha.convert_to(rad));
    
    // Diffusion rate, at the end of the first paragraph of p. 4251
    auto const R_D = (
            species.get_D()[0] * std::pow(sycomore::gamma*G_tau, 2)
        ).convert_to(Hz);
    
    // Equation 5, using definition of R_D
    auto const E1 = [&](int p) { 
        return std::exp(-TR * (R1 + R_D * std::pow(p, 2))); };
    
    // Equation 4, using definition of R_D
    auto const E2 = [&](int p) { 
        return std::exp(-TR * (R2 + R_D*(std::pow(p, 2) + p + 1./3.))); };
    
    // Equation 9
    auto const A = [&](int p) { return 0.5 * (E1(p) - 1.) * (1. + cos_alpha); };
    
    // Equation 10
    auto const B = [&](int p) { return 0.5 * (E1(p) + 1.) * (1. - cos_alpha); };
    
    // Equation 11
    auto const C = [&](int p) { return E1(p) - cos_alpha; };
    
    // Equation 24
    auto const n = [&](int p) { 
        return -E2(-p) * E2(p-1) * std::pow(A(p), 2) * B(p-1) / B(p); };

    // Equation 25
    auto const d = [&](int p) { 
        return A(p)-B(p) + E2(p) * E2(-p-1) * B(p) * C(p+1) / B(p+1); };
    
    // Equation 26 for p=l
    int l = truncation_level;
    auto const e_l = -E2(l) * E2(-l-1) * B(l) * C(l+1) / B(l+1);
    // Equation 23
    auto x_1 = n(truncation_level) / (d(truncation_level) + e_l);
    for(int l=truncation_level-1; l>0; --l)
    {
        x_1 = n(l) / (d(l)+x_1);
    }
    
    // Equation 22
    auto const r_1 = 1/(E2(-1) * B(0)) * x_1 + E2(0) * C(1) / B(1);
    
    // Equation 27, include relaxation term for b_{-1}
    auto const b_0 = 
        -sin_alpha * M0 * (1 - E1(0)) / (A(0) - B(0) + r_1 * E2(-1) * C(0));
    auto const b_minus_1 = -r_1 * E2(-1) * b_0;
    
    // b_0 is the FID mode, b_{-1} is the pre-FID mode.
    return std::make_pair(b_0, -b_minus_1);
}
