// Copyright (c) 2025 Mikko Kuha (University of Jyväskylä).
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef NUCLEUS_GENERATOR_HPP
#define NUCLEUS_GENERATOR_HPP

#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <vector>

#include "nucleon.hpp"
#include "typedefs.hpp"

class nucleus_generator
{
public:
  struct nucleus_params
  {
    // Pb by default
    uint_fast16_t A{208};
    uint_fast16_t ZA{82};
    uint_fast16_t B{208};
    uint_fast16_t ZB{82};
    double RA{6.62435};
    double RB{6.62435};
    double dA{0.5498};
    double dB{0.5498};
    double beta2A{0.0};
    double beta2B{0.0};
    double beta3A{0.0};
    double beta3B{0.0};
    double beta4A{0.0};
    double beta4B{0.0};
    // Overall params
    double min_distance{0.4};
    bool shift_cms{true};
    bool correct_overlap_bias{true};
    // Hotspots
    bool hotspots{false};
    uint_fast16_t n_hotspots{3};
    double hotspot_distr_width{0.5};

    explicit nucleus_params(
        auto A_,
        auto ZA_,
        auto B_,
        auto ZB_,
        auto RA_,
        auto RB_,
        auto dA_,
        auto dB_,
        auto beta2A_,
        auto beta2B_,
        auto beta3A_,
        auto beta3B_,
        auto beta4A_,
        auto beta4B_,
        auto min_distance_,
        auto shift_cms_,
        auto correct_overlap_bias_,
        auto hotspots_,
        auto n_hotspots_,
        auto hotspot_distr_width_) noexcept
        : A(A_),
          ZA(ZA_),
          B(B_),
          ZB(ZB_),
          RA(RA_),
          RB(RB_),
          dA(dA_),
          dB(dB_),
          beta2A(beta2A_),
          beta2B(beta2B_),
          beta3A(beta3A_),
          beta3B(beta3B_),
          beta4A(beta4A_),
          beta4B(beta4B_),
          min_distance(min_distance_),
          shift_cms(shift_cms_),
          correct_overlap_bias(correct_overlap_bias_),
          hotspots(hotspots_),
          n_hotspots(n_hotspots_),
          hotspot_distr_width(hotspot_distr_width_)
    {
    }
  };
  static std::vector<nucleon> generate_nucleus(const nucleus_params &params, const bool &target, const double &mom,
                                               const double &xshift, std::shared_ptr<std::mt19937> random_generator);

protected:
  static coords throw_nucleon_coords(std::shared_ptr<std::mt19937> random_generator, const double &R0, const double &d, const double &beta2, const double &beta3, const double &beta4) noexcept;
  static double throw_radial(const double &random_number) noexcept;
  static void throw_neutrons(std::vector<nucleon> *const nucleus, const uint_fast16_t &Z, std::shared_ptr<std::mt19937> random_generator) noexcept;
  static void throw_hotspots(std::vector<nucleon> *const nucleus, const uint_fast16_t &n_hotspots, const double &hotspot_distr_width, std::shared_ptr<std::mt19937> random_generator) noexcept;
  static bool coords_fit(const coords &co, const std::vector<coords> &other_coords, const double &min_distance) noexcept;

private:
  static std::uniform_real_distribution<double> unif_dist;
  static std::normal_distribution<double> normal_dist;
};

#endif // NUCLEUS_GENERATOR_HPP