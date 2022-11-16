//Copyright (c) 2022 Mikko Kuha

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
        uint_fast16_t A{208};
        uint_fast16_t ZA{82};
        uint_fast16_t B{208};
        uint_fast16_t ZB{82};
         //Overall params
        double min_distance{0.4};
        bool shift_cms{true};
        bool correct_overlap_bias{true};

        explicit nucleus_params
        (
          auto A_, 
          auto ZA_,     
          auto B_,  
          auto ZB_,      
          auto min_distance_,    
          auto shift_cms_,        
          auto correct_overlap_bias_          
        ) noexcept
        : A(A_), 
          ZA(ZA_), 
          B(B_), 
          ZB(ZB_), 
          min_distance(min_distance_), 
          shift_cms(shift_cms_),
          correct_overlap_bias(correct_overlap_bias_)
        {}
    };
    static std::vector<nucleon> generate_nucleus(const nucleus_params & params, const bool &target, const double & mom, 
        const double & xshift, std::shared_ptr<std::mt19937> random_generator, std::shared_ptr<ars> radial_sampler);
protected:
    static coords throw_nucleon_coords(std::shared_ptr<std::mt19937> random_generator, std::shared_ptr<ars> radial_sampler) noexcept;
    static double throw_radial(const double & random_number) noexcept;
    static void throw_neutrons(std::vector<nucleon> *const nucleus, const uint_fast16_t & Z, std::shared_ptr<std::mt19937> random_generator) noexcept;
    static bool coords_fit(const coords& co, const std::vector<coords>& other_coords, const double& min_distance) noexcept;
private:
    static std::uniform_real_distribution<double> unif_dist;
};

#endif // NUCLEUS_GENERATOR_HPP