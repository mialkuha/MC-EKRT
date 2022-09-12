//Copyright (c) 2021 Mikko Kuha

#ifndef NUCLEUS_GENERATOR_HPP
#define NUCLEUS_GENERATOR_HPP

#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <vector>

#include "ars.hpp"
#include "nucleon.hpp"
#include "typedefs.hpp"

class nucleus_generator
{
public:
    struct nucleus_params 
    {
         //Pb by default
        const uint NA{208};
        const uint ZA{82};
        const uint NB{208};
        const uint ZB{82};
         //Overall params
        const spatial min_distance{0.4};
        const bool shift_cms{true};
        const bool correct_overlap_bias{true};
    };
    static std::vector<nucleon> generate_nucleus(const nucleus_params & params, const bool &target, const momentum & mom, 
        const spatial & xshift, std::shared_ptr<std::mt19937> random_generator, std::shared_ptr<ars> radial_sampler);
protected:
    static coords&& throw_nucleon_coords(std::shared_ptr<std::mt19937> random_generator, std::shared_ptr<ars> radial_sampler) noexcept;
    static double throw_radial(const double & random_number) noexcept;
    static void throw_neutrons(std::vector<nucleon> *const nucleus, const uint & Z, std::shared_ptr<std::mt19937> random_generator) noexcept;
    static bool coords_fit(const coords& co, const std::vector<coords>& other_coords, const spatial& min_distance) noexcept;
private:
    static std::uniform_real_distribution<double> unif_dist;
};

#endif // NUCLEUS_GENERATOR_HPP