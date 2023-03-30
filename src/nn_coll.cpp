//Copyright (c) 2021 Mikko Kuha

#include "nn_coll.hpp"

void nn_coll::reduce_energy() noexcept
{
    //NO boosts required when massles nucleons and collinear partons
    auto x1 = 0.0, x2 = 0.0;
    for (auto s : this->dijets)
    {
        x1 += (s.kt / this->sqrt_s) * (exp(s.y1) + exp(s.y2));
        x2 += (s.kt / this->sqrt_s) * (exp(-s.y1) + exp(-s.y2));
    }
    

    this->projectile->mom *= 1.0 - x1;
    this->target->mom *= 1.0 - x2;
}

void nn_coll::wound() noexcept
{
    this->projectile->wounded = true;
    this->target->wounded = true;
}

void nn_coll::calculate_xsects(const double &sigma_jet, const std::function<double(const double&)> &Tpp, const double &b2) noexcept
{
    double max_tot=1., max_inel=1., eff_tot=1., eff_inel=1.;

    max_tot = 2. * M_PI * (1. - exp(-0.5 * sigma_jet * Tpp(0.)));
    eff_tot = 2. * M_PI * (1. - exp(-0.5 * sigma_jet * Tpp(b2)));
    max_inel = M_PI * (1. - exp(-sigma_jet * Tpp(0.)));
    eff_inel = M_PI * (1. - exp(-sigma_jet * Tpp(b2)));

    this->max_inel_xsect = max_inel; //sigma_inel
    this->max_tot_xsect = max_tot; //sigma_tot
    this->effective_inel_xsect = eff_inel; //sigma_inel
    this->effective_tot_xsect = eff_tot; //sigma_tot
}

void nn_coll::push_end_states_to_collider_frame() noexcept
{
    //Calculate transformation btw events CMS and collider frame
    //double rap = atanh( (this->projectile->mom - this->target->mom) 
    //                    / (this->projectile->mom + this->target->mom));
    double rap = 0.5*log(this->projectile->mom / this->target->mom);

    //Push all end states with the transformation
    for (auto &s : this->dijets)
    {
        s.y1 += rap;
        s.y2 += rap;
    }
}

void nn_coll::reduce_energy_and_push_end_states_to_collider_frame() noexcept
{
    //Calculate transformation btw events CMS and collider frame
    //double rap = atanh( (this->projectile->mom - this->target->mom) 
    //                    / (this->projectile->mom + this->target->mom));
    double rap = 0.5*log(this->projectile->mom / this->target->mom);
    auto x1 = 0.0, x2 = 0.0;
    for (auto s : this->dijets)
    {
        x1 += (s.kt / this->sqrt_s) * (exp(s.y1) + exp(s.y2));
        x2 += (s.kt / this->sqrt_s) * (exp(-s.y1) + exp(-s.y2));
    }

    this->projectile->mom *= 1.0 - x1;
    this->target->mom *= 1.0 - x2;

    //Push all end states with the transformation
    for (auto &s : this->dijets)
    {
        s.y1 += rap; // y1_LAB = y1_CMS + rap
        s.y2 += rap;
    }
}