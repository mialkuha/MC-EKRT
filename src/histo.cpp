//Copyright (c) 2022 Mikko Kuha

#include "histo.hpp"

void histo_1d::add
(
    const double &y
) noexcept
{
    if (y <= this->xs.front())
    {
        this->underf++;
        this->total_counts++;
        return;
    }
    
    uint16_t y_index = 0;

    for 
    (
        auto x_it = ++this->xs.begin();
        x_it!=this->xs.end(); 
        x_it++
    )
    {
        if (y <= *x_it)
        {
            this->counts[y_index]++;
            this->total_counts++;
            return;
        }
        y_index++;
    }

    this->overf++;
    this->total_counts++;
}