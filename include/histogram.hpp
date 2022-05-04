//Copyright (c) 2022 Mikko Kuha

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>

namespace histogram
{
    class histo_1d
    {
    public:
        histo_1d
        (
            std::vector<double> &xs_
        ) noexcept : 
            xs(xs_),
            counts(xs_.size(), 0.0)
        {}
        histo_1d
        (
            std::vector<double> &xs_,
            std::vector<double> &ys_
        ) noexcept : 
            xs(xs_),
            counts(xs_.size()-1, 0.0)
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
            const histogram::histo_1d &other
        ) noexcept -> void;

        auto get_histo() noexcept;

    protected:
    private:
        const std::vector<double> xs;
        const std::vector<double> counts;
        uint64_t total_counts;
        double underf{0.0};
        double overf{0.0};
    };
};

#endif // HISTOGRAM_H
