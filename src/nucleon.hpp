//Copyright (c) 2021 Mikko Kuha

#ifndef NUCLEON_HPP
#define NUCLEON_HPP

#include <cmath>
#include <iostream>

#include "typedefs.hpp"

class nucleon
{
public:
    nucleon(coords co_, const double &mom_) noexcept 
        : co(std::move(co_)), mom(mom_), wounded(false), is_neutron(false) {};

    nucleon(coords co_, const double & mom_, const bool & wounded_, const bool & is_neutron_) noexcept 
        : co(std::move(co_)), mom(mom_), wounded(wounded_), is_neutron(is_neutron_) {};

    //nucleon& operator=(const nucleon &) = default;
    //nucleon(const nucleon &) = default;

    double calculate_bsquared(const nucleon &other) noexcept
    {
        return std::pow(other.co.x - this->co.x, 2) + std::pow(other.co.y - this->co.y, 2);
    }

    coords co{0,0,0};
    double mom{0};
    bool wounded{false};
    bool is_neutron{false};
protected:
private:
};

#endif // NUCLEON_HPP
