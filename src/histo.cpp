//Copyright (c) 2022 Mikko Kuha

#include "histo.hpp"

auto histo_1d::add
(
    const double &y
) noexcept -> void
{
    if (y <= this->xs[0])
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

auto histo_1d::add
(
    const std::vector<double> &news
) noexcept -> void
{
    uint64_t news_tot = news.size();

    if (news_tot <= 0)
    {
        return;
    }
    else if (news_tot == 1)
    {
        this->add(news[0]);
        return;
    }

    uint64_t news_index = 0;
    std::vector<double> to_add{news};
    std::sort(to_add.begin(), to_add.end());

    while (to_add[news_index] <= this->xs[0])
    {
        this->underf++;
        this->total_counts++;
        news_index++;
        if (news_index == news_tot)
        {
            return;
        }
    }
    
    uint16_t counts_index = 0;

    for 
    (
        auto x_it = ++this->xs.begin();
        x_it!=this->xs.end(); 
        x_it++
    )
    {
        while (to_add[news_index] <= *x_it)
        {
            this->counts[counts_index]++;
            this->total_counts++;
            news_index++;
            if (news_index == news_tot)
            {
                return;
            }
        }
        counts_index++;
    }

    while (news_index < news_tot)
    {
        this->overf++;
        this->total_counts++;
        news_index++;
    }
}

auto histo_1d::add
(
    const histo_1d &other
) noexcept -> void
{
    if (!(*this == other))
    {
        std::cout<<"Incompatible histograms!"<<std::endl;
        return;
    }

    this->underf += other.underf;
    this->overf += other.overf;
    for (uint16_t i = 0; i<this->counts.size(); i++)
    {
        this->counts[i] += other.counts[i];
    }
    this->total_counts += other.total_counts;
}

auto histo_1d::get_histo() const noexcept
{   
    std::vector<double> ret_histo{this->counts};
    auto total = this->total_counts;
    for (auto & r :ret_histo)
    {
        r /= total;
    }

    return std::make_tuple
           (
               std::move(ret_histo),
               this->underf / total,
               this->overf / total,
               this->xs, 
               this->total_counts
           );
}

auto operator==
(
    const histo_1d& c1, 
    const histo_1d& c2
) noexcept -> bool
{
    bool ret = (c1.xs.size() == c2.xs.size());
    ret *= (c1.counts.size() == c2.counts.size());
    for (uint16_t i = 0; i<c1.xs.size(); i++)
    {
        ret *= (c1.xs[i] == c2.xs[i]);
    }
    return ret;
}

int main()
{
    std::vector<double> xs{0.0,0.5,1.0,1.5,2.5};
    histo_1d h1{xs};
    auto [ ys, uf, of, xx, tot ] = h1.get_histo();
    
    for (auto x : ys) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto x : xx) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf<<std::endl;
    std::cout<<of<<std::endl;
    std::cout<<tot<<std::endl;

    h1.add(0.1);
    h1.add(std::vector<double>({1.1,0.2,-2}));
    auto [ ys1, uf1, of1, xx1, tot1 ] = h1.get_histo();
    
    for (auto x : ys1) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto x : xx1) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf1<<std::endl;
    std::cout<<of1<<std::endl;
    std::cout<<tot1<<std::endl;

    return 0;
}