//Copyright (c) 2021 Mikko Kuha

#include "nucleon.hpp"

nucleon::nucleon(coords co_, const momentum &mom_) noexcept : co(std::move(co_)), mom(mom_), wounded(false), is_neutron(false)
{
}

spatial nucleon::calculate_bsquared(const nucleon &other) noexcept
{
    return std::pow(other.co.x - this->co.x, 2) + std::pow(other.co.y - this->co.y, 2);
}
