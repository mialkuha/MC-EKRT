// Copyright (c) 2025 Mikko Kuha.
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software 
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef LINEAR_INTERPOLATOR_H
#define LINEAR_INTERPOLATOR_H

#include <iostream>
#include <vector>

class linear_interpolator
{
public:
    linear_interpolator(std::vector<double> xs_, std::vector<double> ys_) noexcept : xs(std::move(xs_)), ys(std::move(ys_)) {}
    double value_at(const double &x) const noexcept;
            
    linear_interpolator& operator=(const linear_interpolator& rhs)
    {
        this->xs = rhs.xs;
        this->ys = rhs.ys;
        return *this;
    }

protected:
private:
    std::vector<double> xs;
    std::vector<double> ys;
};

#endif // LINEAR_INTERPOLATOR_H
