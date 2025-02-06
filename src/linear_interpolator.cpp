// Copyright (c) 2025 Mikko Kuha.
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software 
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#include "linear_interpolator.hpp"

double linear_interpolator::value_at(const double &x) const noexcept
{
    if (x < this->xs.front() || x > this->xs.back())
    {
        std::cout << "LINEAR INTERPOLATOR OUT OF BOUNDS AT x="<<x<< std::endl;
        return 0.0;
    }

    if (x == this->xs.front())
    {
        return this->ys.front();
    }

    size_t nearest_upper = 0;

    for (size_t i = 0; i < this->xs.size(); i++)
    {
        if (x <= this->xs[i])
        {
            nearest_upper = i;
            break;
        }
    }

    if (x == this->xs[nearest_upper])
    {
        return this->ys[nearest_upper];
    }
    return (this->ys[nearest_upper - 1] * (this->xs[nearest_upper] - x) + this->ys[nearest_upper] * (x - this->xs[nearest_upper - 1])) / (this->xs[nearest_upper] - this->xs[nearest_upper - 1]);
}
