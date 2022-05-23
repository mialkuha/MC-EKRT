//Copyright (c) 2021 Mikko Kuha

#ifndef NUCLEON_HPP
#define NUCLEON_HPP

#include <cmath>
#include <iostream>

#include "typedefs.hpp"

class nucleon
{
public:
    nucleon(coords co_, const momentum & mom_) noexcept;
    nucleon(coords co_, const momentum & mom_, const bool & wounded_, const bool & is_neutron_)
        noexcept : co(std::move(co_)), mom(mom_), wounded(wounded_), is_neutron(is_neutron_) {};

    spatial calculate_bsquared(const nucleon & other) noexcept;

    coords co{0,0,0};
    momentum mom{0};
    bool wounded{false};
    bool is_neutron{false};
protected:
private:
};

#endif // NUCLEON_HPP
