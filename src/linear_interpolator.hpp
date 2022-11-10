//Copyright (c) 2021 Mikko Kuha

#ifndef LINEAR_INTERPOLATOR_H
#define LINEAR_INTERPOLATOR_H

#include <iostream>
#include <vector>

class linear_interpolator
{
public:
    linear_interpolator(std::vector<double> xs_, std::vector<double> ys_) noexcept : xs(std::move(xs_)), ys(std::move(ys_)) {}
    double value_at(const double &x) noexcept;

protected:
private:
    const std::vector<double> xs;
    const std::vector<double> ys;
};

#endif // LINEAR_INTERPOLATOR_H
