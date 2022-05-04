//Copyright (c) 2022 Mikko Kuha

#include "histogram.hpp"

auto histogram::histo_1d::add
(
    const double &y
) noexcept -> void
{
    if (y <= this->xs.front())
    {
        this->underf++;
        total_counts++;
        return;
    }

    uint16_t xs_size = this->xs.size();

    for 
    (
        auto x_it = ++this->xs.begin(),
             y_it = this->counts.begin();
        x_it!=this->xs.end(); 
        x_it++, y_it++
    )
    {
        if (y <= *x_it)
        {
            (*y_it)++;
            total_counts++;
            return;
        }
    }

    this->overf++;
    total_counts++;
}