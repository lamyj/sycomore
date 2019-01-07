#ifndef _281032dc_2a2d_4545_8430_64b745cb4f19
#define _281032dc_2a2d_4545_8430_64b745cb4f19

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sycomore/Grid.h"
#include "sycomore/magnetization.h"
#include "sycomore/Pulse.h"
#include "sycomore/Species.h"
#include "sycomore/sycomore.h"
#include "sycomore/TimeInterval.h"

namespace sycomore
{

class Model
{
public:
    /// @brief Relative B1
    Real B1 = 1.;

    /**
     * @brief Create a new model from an initial state (defined by a species and
     * its initial real magnetizations) and a set of named time intervals.
     * magnetization.
     */
    Model(
        Species const & species, Magnetization const & magnetization,
        std::vector<std::pair<std::string, TimeInterval>> const & time_intervals);

    /// @brief Return the dimension index map.
    std::map<std::string, size_t> const & dimensions() const;

    /// @brief Return the models for all species
    std::map<std::string, TimeInterval> const & time_intervals() const;

    /// @brief Return the threshold magnetization for clean-up (default to 0).
    Real epsilon() const;

    /// @brief Set the threshold magnetization for clean-up.
    void epsilon(Real const & value);

    /// @brief Apply an RF pulse to the model, modulated by the relative B1.
    void apply_pulse(Pulse const & pulse);

    /// @brief Apply a time interval to the model.
    void apply_time_interval(std::string const & name);

    /// @brief Return the complex magnetizations.
    Grid const & grid() const;

    /**
     * @brief Return the isochromat for the given configurations.
     *
     * If no configuration is specified, all configurations are used.
     */
    Magnetization isochromat(std::vector<Index> const & configurations={}) const;

private:
    /// @brief Species of the model.
    Species _species;

    /// @brief Dimension index of time interval
    std::map<std::string, size_t> _dimensions;

    /// @brief Named time intervals.
    std::map<std::string, TimeInterval> _time_intervals;

    /// @brief Precision to cull low-populated states
    Real _epsilon_squared;

    /// @brief the configuration models, one for each species.
    Grid _grid;

    /// @brief Bounding box of the occupied configuration for each species.
    std::pair<Index, Shape> _bounding_box;

    void _cleanup(Grid & model, std::pair<Index, Shape> const & bounding_box);
};

}

#endif // _281032dc_2a2d_4545_8430_64b745cb4f19
