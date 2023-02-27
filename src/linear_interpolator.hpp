//Copyright (c) 2023 Mikko Kuha

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
