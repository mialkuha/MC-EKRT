//Copyright (c) 2022 Mikko Kuha

#ifndef HISTO_H
#define HISTO_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>

class histo_1d
{
public:
    histo_1d
    (
        const std::vector<double> &xs_
    ) noexcept : 
        xs(xs_),
        counts(xs_.size(), 0.0),
        total_counts(0),
        underf(0),
        overf(0)
    {}
    histo_1d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_
    ) noexcept : 
        xs(xs_),
        counts(xs_.size()-1, 0.0),
        total_counts(0),
        underf(0),
        overf(0)
    { this->add(ys_); }

    auto add
    (
        const double &y
    ) noexcept -> void;

    auto add
    (
        const std::vector<double> &ys
    ) noexcept -> void;

    auto add
    (
        const histo_1d &other
    ) noexcept -> void;

    auto get_histo() const noexcept;

    friend bool operator==
    (
        const histo_1d& c1, 
        const histo_1d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    std::vector<double> counts;
    uint64_t total_counts{0};
    double underf{0.0};
    double overf{0.0};
};

#endif // HISTO_H
