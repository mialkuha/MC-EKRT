//Copyright (c) 2022 Mikko Kuha

#ifndef HISTO_H
#define HISTO_H

#include <cstdint>
#include <vector>

class histo_1d
{
public:
    histo_1d
    (
        const std::vector<double> &xs_
    ) noexcept : 
        xs(xs_),
        counts(xs_.size(), 0.0)
    {}
    histo_1d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_
    ) noexcept : 
        xs(xs_),
        counts(xs_.size()-1, 0.0)
    { this->add(ys_); }

    void add
    (
        const double &y
    ) noexcept;

    auto add
    (
        const std::vector<double> &ys
    ) noexcept -> void;

    auto add
    (
        const histo_1d &other
    ) noexcept -> void;

    auto get_histo() const noexcept;

protected:
private:
    const std::vector<double> xs;
    std::vector<double> counts;
    uint64_t total_counts;
    double underf{0.0};
    double overf{0.0};
};

#endif // HISTO_H
