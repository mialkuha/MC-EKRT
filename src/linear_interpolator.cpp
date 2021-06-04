//Copyright (c) 2021 Mikko Kuha

#include "linear_interpolator.hpp"

double linear_interpolator::value_at(const double &x) noexcept
{
    if (x < this->xs.front() || x > this->xs.back())
    {
        std::cout << "OUT OF BOUNDS" << std::endl;
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
