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
        counts(xs_.size()-1, 0.0),
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

class histo_2d
{
public:
    histo_2d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_
    ) noexcept : 
        xs(xs_),
        ys(ys_),
        counts
        (
            xs_.size()-1, 
            std::vector<double>
            (
                ys_.size()-1, 
                0.0
            )
        ),
        total_counts(0),
        underf({0.0, 0.0}),
        overf({0.0, 0.0})
    {}

    auto add
    (
        const std::tuple<double, double> &x
    ) noexcept -> void;

    auto add
    (
        std::vector<std::tuple<double, double> > news
    ) noexcept -> void;

    auto add
    (
        const histo_2d &other
    ) noexcept -> void;

    auto get_histo() const noexcept;

    auto project_1d(const bool project_ys) const noexcept;

    friend bool operator==
    (
        const histo_2d& c1, 
        const histo_2d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    const std::vector<double> ys;
    std::vector<std::vector<double> > counts;
    uint64_t total_counts{0};
    std::tuple<double, double> underf{0.0, 0.0};
    std::tuple<double, double> overf{0.0, 0.0};
};

class histo_3d
{
public:
    histo_3d
    (
        const std::vector<double> &xs_,
        const std::vector<double> &ys_,
        const std::vector<double> &zs_
    ) noexcept : 
        xs(xs_),
        ys(ys_),
        zs(zs_),
        counts
        (
            xs_.size()-1, 
            std::vector<std::vector<double> >
            (
                ys_.size()-1, 
                std::vector<double>
                (
                    zs_.size()-1, 
                    0.0
                )
            )
        ),
        total_counts(0),
        underf({0.0, 0.0, 0.0}),
        overf({0.0, 0.0, 0.0})
    {}

    auto add
    (
        const std::tuple<double, double, double> &x
    ) noexcept -> void;

    auto add
    (
        std::vector<std::tuple<double, double, double> > news
    ) noexcept -> void;

    auto add
    (
        const histo_3d &other
    ) noexcept -> void;

    auto get_histo() const noexcept;

    auto project_1d(const uint8_t dim_left) const noexcept;

    friend bool operator==
    (
        const histo_3d& c1, 
        const histo_3d& c2
    ) noexcept;

protected:
private:
    const std::vector<double> xs;
    const std::vector<double> ys;
    const std::vector<double> zs;
    std::vector<std::vector<std::vector<double> > > counts;
    uint64_t total_counts{0};
    std::tuple<double, double, double> underf{0.0, 0.0, 0.0};
    std::tuple<double, double, double> overf{0.0, 0.0, 0.0};
};

#endif // HISTO_H
