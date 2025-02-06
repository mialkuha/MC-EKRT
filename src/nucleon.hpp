// Copyright (c) 2025 Mikko Kuha.
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software 
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef NUCLEON_HPP
#define NUCLEON_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "typedefs.hpp"

class nucleon
{
public:
    nucleon(coords co_, const double &mom_, const int_fast16_t &index_) noexcept 
        : co(std::move(co_)), mom(mom_), wounded(false), is_neutron(false), index(index_) {};

    nucleon(coords co_, const double & mom_, const bool & wounded_, const bool & is_neutron_, const int_fast16_t &index_) noexcept 
        : co(std::move(co_)), mom(mom_), wounded(wounded_), is_neutron(is_neutron_), index(index_) {};

    //nucleon& operator=(const nucleon &) = default;
    //nucleon(const nucleon &) = default;

    double calculate_bsquared(const nucleon &other) noexcept
    {
        return std::pow(other.co.x - this->co.x, 2) + std::pow(other.co.y - this->co.y, 2);
    }

    double calculate_bsquared_hotspot(const nucleon &other) noexcept
    {
        double bsquared = 100.0;
        for (auto pro : this->hotspots)
        {
            for (auto tar : other.hotspots)
            {
                double dummy = std::pow(pro.co.x - tar.co.x, 2) + std::pow(pro.co.y - tar.co.y, 2);
                if (dummy < bsquared)
                {
                    bsquared = dummy;
                }
            }
        }
        return bsquared;
    }

    coords co{0,0,0};
    double mom{0};
    std::vector<hotspot_info> hotspots;
    bool wounded{false};
    bool is_neutron{false};
    int_fast16_t index{0};

    friend bool operator==(const nucleon &lhs, const nucleon &rhs)
    {
        return (lhs.co == rhs.co)&&
            (lhs.mom == rhs.mom)&&
            (lhs.wounded == rhs.wounded)&&
            (lhs.is_neutron == rhs.is_neutron)&&
            (lhs.index == rhs.index);
    }
    
protected:
private:
};

#endif // NUCLEON_HPP
