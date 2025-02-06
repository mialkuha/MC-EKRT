// Copyright (c) 2025 Mikko Kuha (University of Jyväskylä).
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#include "nucleus_generator.hpp"

std::uniform_real_distribution<double> nucleus_generator::unif_dist = std::uniform_real_distribution<double>(0.0, 1.0);
std::normal_distribution<double> nucleus_generator::normal_dist = std::normal_distribution<double>(0, 0);

std::vector<nucleon> nucleus_generator::generate_nucleus(const nucleus_params &params, const bool &target, const double &mom,
                                                         const double &xshift, std::shared_ptr<std::mt19937> random_generator)
{
    std::vector<nucleon> generated_nucleus;
    std::vector<coords> generated_coords;
    auto N = (target) ? params.B : params.A;
    auto Z = (target) ? params.ZB : params.ZA;
    auto R0 = (target) ? params.RB : params.RA;
    auto d = (target) ? params.dB : params.dA;
    auto beta2 = (target) ? params.beta2B : params.beta2A;
    auto beta3 = (target) ? params.beta3B : params.beta3A;
    auto beta4 = (target) ? params.beta4B : params.beta4A;
    generated_nucleus.reserve(N);
    generated_coords.reserve(N);

    if (N == 1)
    {
        generated_nucleus.emplace_back(coords({xshift, 0.0, 0.0}), mom, 0);
        if (params.hotspots)
        {
            nucleus_generator::throw_hotspots(&generated_nucleus, params.n_hotspots, params.hotspot_distr_width, random_generator);
        }
        return generated_nucleus;
    }

    coords new_coords;
    bool coords_do_fit;
    do // while (generated_coords.size() < N);
    {
        do // while (!coords_do_fit);
        {
            try
            {
                new_coords = nucleus_generator::throw_nucleon_coords(random_generator, R0, d, beta2, beta3, beta4);
                coords_do_fit = nucleus_generator::coords_fit(new_coords, generated_coords, params.min_distance);
                // std::cout<<"coords "<<new_coords.x<<' '<<new_coords.y<<' '<<new_coords.z<<" fit? "<<coords_do_fit<<std::endl;
                if (params.correct_overlap_bias && !coords_do_fit)
                {
                    break;
                }
            }
            catch (const std::exception &e)
            {
                std::cout << e.what() << " in nucleus_generator, trying again" << std::endl;
                coords_do_fit = false;
            }
        } while (!coords_do_fit);

        if (params.correct_overlap_bias && !coords_do_fit)
        {
            generated_coords.clear();
        }
        generated_coords.push_back(std::move(new_coords));

    } while (generated_coords.size() < N);

    coords com{0.0, 0.0, 0.0}; // center of mass

    if (params.shift_cms)
    {
        for (const auto &co : generated_coords)
        {
            com += co;
        }
        com /= static_cast<double>(generated_coords.size());
    }

    int_fast16_t index = 0;
    for (auto &co : generated_coords)
    {
        if (params.shift_cms)
        {
            co -= com; // Shift all nucleons so that CoM is (0,0,0)
        }
        co.x += xshift;
        generated_nucleus.emplace_back(co, mom, index++);
    }

    nucleus_generator::throw_neutrons(&generated_nucleus, Z, random_generator);
    if (params.hotspots)
    {
        nucleus_generator::throw_hotspots(&generated_nucleus, params.n_hotspots, params.hotspot_distr_width, random_generator);
    }

    return generated_nucleus;
}

coords nucleus_generator::throw_nucleon_coords(std::shared_ptr<std::mt19937> random_generator, const double &R0, const double &d, const double &beta2, const double &beta3, const double &beta4) noexcept
{
    double ws = 0.0, rand4 = 0.0, x = 0.0, y = 0.0, z = 0.0;
    do
    {
        double rand1 = nucleus_generator::unif_dist(*random_generator);
        double rand2 = nucleus_generator::unif_dist(*random_generator);
        double rand3 = nucleus_generator::unif_dist(*random_generator);
        rand4 = nucleus_generator::unif_dist(*random_generator);

        x = (3.0 * R0) * rand1 - (3.0 * R0 / 2.0);
        y = (3.0 * R0) * rand2 - (3.0 * R0 / 2.0);
        z = (3.0 * R0) * rand3 - (3.0 * R0 / 2.0);

        double r, cstheta, cstheta2, R, Rmax;

        r = std::sqrt(x * x + y * y + z * z);

        if (r > 0)
        {
            cstheta2 = z * z / (r * r);
            cstheta = z / r;
        }
        else
        {
            cstheta2 = 0.0;
            cstheta = 0.0;
        }
        R = R0 * (1 + beta2 * 1.0 / 4.0 * std::sqrt(5.0 / M_PI) * (3.0 * cstheta2 - 1.0) + beta3 * 0.25 * std::sqrt(7.0 / M_PI) * (5.0 * cstheta2 * cstheta - 3.0 * cstheta) + beta4 * 3.0 / 16.0 * std::sqrt(1.0 / M_PI) * (35.0 * cstheta2 * cstheta2 - 30.0 * cstheta2 + 3.0));

        Rmax = R0 * (1 + std::abs(beta2) * 1.0 / 2.0 * std::sqrt(5.0 / M_PI) + std::abs(beta3) * 0.5 * std::sqrt(7.0 / M_PI) + std::abs(beta4) * 3.0 / 2.0 * std::sqrt(1.0 / M_PI));

        ws = 1.0 / (1.0 + std::exp((r - R) / d)) * (1.0 + std::exp((-Rmax) / d));

    } while (rand4 > ws);

    return coords({x, y, z});
}

void nucleus_generator::throw_neutrons(std::vector<nucleon> *const nucleus, const uint_fast16_t &Z, std::shared_ptr<std::mt19937> random_generator) noexcept
{
    for (uint_fast16_t i = 0, iz = 0; i < nucleus->size(); ++i)
    {
        double frac = static_cast<double>(Z - iz) / static_cast<double>(nucleus->size() - i);
        double rn = nucleus_generator::unif_dist(*random_generator);
        if (rn < frac)
        {
            nucleus->at(i).is_neutron = false;
            ++iz;
        }
        else
        {
            nucleus->at(i).is_neutron = true;
        }
    }
}

void nucleus_generator::throw_hotspots(std::vector<nucleon> *const nucleus, const uint_fast16_t &n_hotspots, const double &hotspot_distr_width, std::shared_ptr<std::mt19937> random_generator) noexcept
{
    for (uint_fast16_t i = 0; i < nucleus->size(); ++i)
    {
        nucleus->at(i).hotspots = std::vector<hotspot_info>();
        for (uint_fast16_t j = 0; j < n_hotspots; ++j)
        {
            auto param = std::normal_distribution<double>::param_type{0., hotspot_distr_width};
            coords ds = {normal_dist(*random_generator, param), normal_dist(*random_generator, param), 0.0};
            auto hs_coords = nucleus->at(i).co + ds;

            nucleus->at(i).hotspots.push_back(hotspot_info{hs_coords});
        }
    }
}

bool nucleus_generator::coords_fit(const coords &co, const std::vector<coords> &other_coords, const double &min_distance) noexcept
{
    if (min_distance <= 0)
    {
        return true;
    }

    const double md2 = pow(min_distance, 2);

    for (const auto &oco : other_coords)
    {
        if ((co - oco).mag2() < md2)
        {
            // std::cout<<co.x<<' '<<co.y<<' '<<co.z<<" is too near to "<<oco.x<<' '<<oco.y<<' '<<oco.z<<'\n';
            return false;
        }
    }
    return true;
}
