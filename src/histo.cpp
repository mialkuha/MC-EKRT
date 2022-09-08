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
    const std::tuple<double, double> &y_
) noexcept -> void
{
    auto [ y , w ] = y_;

    if (y <= this->xs[0])
    {
        this->underf += w;
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
            this->counts[y_index] += w;
            this->total_counts++;
            return;
        }
        y_index++;
    }

    this->overf += w;
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
    const std::vector<std::tuple<double, double> > &news
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
    std::vector<std::tuple<double, double> > to_add{news};
    std::sort
    (
        to_add.begin(), 
        to_add.end(),
        [](std::tuple<double, double> a, std::tuple<double, double> b)
            {
                auto [ a0, ph1 ] = a;
                auto [ b0, ph2 ] = b;
                return a0 < b0;
            }
    );

    while (std::get<0>(to_add[news_index]) <= this->xs[0])
    {
        this->underf += std::get<1>(to_add[news_index]);
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
        while (std::get<0>(to_add[news_index]) <= *x_it)
        {
            this->counts[counts_index] += std::get<1>(to_add[news_index]);
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
        this->overf += std::get<1>(to_add[news_index]);
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

auto histo_2d::add
(
    const std::tuple<double, double> &x
) noexcept -> void
{
    auto [ new_x, new_y ] = x;

    {//Overflow or underflow ?
        if (new_x <= this->xs[0])
        {
            std::get<0>(this->underf)++;
            this->total_counts++;
            return;
        }
        else if (new_y <= this->ys[0])
        {
            std::get<1>(this->underf)++;
            this->total_counts++;
            return;
        }
        else if (new_x > this->xs.back())
        {
            std::get<0>(this->overf)++;
            this->total_counts++;
            return;
        }
        else if (new_y > this->ys.back())
        {
            std::get<1>(this->overf)++;
            this->total_counts++;
            return;
        }
    }
    
    uint16_t x_index = 0, y_index = 0;
    for 
    (
        auto x_it = ++this->xs.begin();
        x_it!=this->xs.end(); 
        x_it++
    )
    {
        if (new_x <= *x_it)
        {
            for 
            (
                auto y_it = ++this->ys.begin();
                y_it!=this->ys.end(); 
                y_it++
            )
            {
                if (new_y <= *y_it)
                {
                    this->counts[x_index][y_index]++;
                    this->total_counts++;
                    return;
                }
                y_index++;
            }
            std::cout<<"ERROR, y-bin not found"<<std::endl;
        }
        x_index++;
    }
    std::cout<<"ERROR, x-bin not found"<<std::endl;
}

auto histo_2d::add
(
    std::vector<std::tuple<double, double> > news
) noexcept -> void
{
    //Remove all overflow and underflow
    auto uoflow = [&](std::tuple<double, double> n)
        {
            if (std::get<0>(n) <= this->xs[0])
            {
                std::get<0>(this->underf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<1>(n) <= this->ys[0])
            {
                std::get<1>(this->underf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<0>(n) > this->xs.back())
            {
                std::get<0>(this->overf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<1>(n) > this->ys.back())
            {
                std::get<1>(this->overf)++;
                this->total_counts++;
                return true;
            }
            else
            {
                return false;
            }
        }; 
    std::erase_if(news, uoflow);
    //news.erase
    //(  
    //    std::remove_if(news.begin(), news.end(), uoflow), 
    //    news.end()
    //);

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
    std::vector<std::tuple<double, double> > to_add{news};
    std::vector<double> to_add_1bin{};
    to_add_1bin.reserve(to_add.size());

    //Sort news by x coordinate
    std::sort
    (
        to_add.begin(), 
        to_add.end(),
        [](std::tuple<double, double> a, std::tuple<double, double> b)
            {
                auto [ a0, ph1 ] = a;
                auto [ b0, ph2 ] = b;
                return a0 < b0;
            }
    );
    
    uint16_t x_index = 0, y_index = 0;
    uint64_t x_bin_tot = 0, to_add_index = 0;
    
    //Go through all the x-bins
    for 
    (
        auto x_it = ++this->xs.begin();
        x_it!=this->xs.end(); 
        x_it++
    )
    {
        //Add all the coords that belong in one x-bin at the same time
        while (std::get<0>(to_add[news_index]) <= *x_it)
        {
            //Stage all coordinates that are in this x-bin
            to_add_1bin.push_back(std::get<1>(to_add[news_index]));
            news_index++;
            if (news_index == news_tot)
            {
                break;
            }
        }
        if (to_add_1bin.size()>0)
        {
            to_add_index = 0;
            x_bin_tot = to_add_1bin.size();

            //Sort the coords by y coordinate
            std::sort(to_add_1bin.begin(), to_add_1bin.end());
            y_index = 0;

            //Go through all the y-bins
            for 
            (
                auto y_it = ++this->ys.begin();
                y_it!=this->ys.end(); 
                y_it++
            )
            {
                while (to_add_1bin[to_add_index] <= *y_it)
                {
                    this->counts[x_index][y_index]++;
                    this->total_counts++;
                    to_add_index++;
                    if (to_add_index == x_bin_tot)
                    {
                        goto added_staged;
                    }
                }
                y_index++;
            }
            //All the y-bins checked
            std::cout<<"ERROR, y-bin not found"<<std::endl;
        }
    added_staged:        
        if (news_index == news_tot)
        {
            //All the new coords are added
            return;
        }
        to_add_1bin.clear();
        x_index++;
    }
    //All the x-bins checked
    std::cout<<"ERROR, x-bin not found"<<std::endl;
}

auto histo_2d::add
(
    const histo_2d &other
) noexcept -> void
{
    if (!(*this == other))
    {
        std::cout<<"Incompatible histograms!"<<std::endl;
        return;
    }

    std::get<0>(this->underf) += std::get<0>(other.underf);
    std::get<1>(this->underf) += std::get<1>(other.underf);

    std::get<0>(this->overf) += std::get<0>(other.overf);
    std::get<1>(this->overf) += std::get<1>(other.overf);

    for (uint16_t i = 0; i<this->counts.size(); i++)
    {
        for (uint16_t j = 0; j<this->counts[0].size(); j++)
        {
            this->counts[i][j] += other.counts[i][j];
        }
    }

    this->total_counts += other.total_counts;
}

auto histo_3d::add
(
    const std::tuple<double, double, double> &x
) noexcept -> void
{
    auto [ new_x, new_y, new_z ] = x;

    {//Overflow or underflow ?
        if (new_x <= this->xs[0])
        {
            std::get<0>(this->underf)++;
            this->total_counts++;
            return;
        }
        else if (new_y <= this->ys[0])
        {
            std::get<1>(this->underf)++;
            this->total_counts++;
            return;
        }
        else if (new_z <= this->zs[0])
        {
            std::get<2>(this->underf)++;
            this->total_counts++;
            return;
        }
        else if (new_x > this->xs.back())
        {
            std::get<0>(this->overf)++;
            this->total_counts++;
            return;
        }
        else if (new_y > this->ys.back())
        {
            std::get<1>(this->overf)++;
            this->total_counts++;
            return;
        }
        else if (new_z > this->zs.back())
        {
            std::get<2>(this->overf)++;
            this->total_counts++;
            return;
        }
    }
    
    uint16_t x_index = 0, y_index = 0, z_index = 0;
    for 
    (
        auto x_it = ++this->xs.begin();
        x_it!=this->xs.end(); 
        x_it++
    )
    {
        if (new_x <= *x_it)
        {
            for 
            (
                auto y_it = ++this->ys.begin();
                y_it!=this->ys.end(); 
                y_it++
            )
            {
                if (new_y <= *y_it)
                {
                    for 
                    (
                        auto z_it = ++this->zs.begin();
                        z_it!=this->zs.end(); 
                        z_it++
                    )
                    {
                        if (new_z <= *z_it)
                        {
                            this->counts[x_index][y_index][z_index]++;
                            this->total_counts++;
                            return;
                        }
                        z_index++;
                    }
                    std::cout<<"ERROR, z-bin not found"<<std::endl;
                }
                y_index++;
            }
            std::cout<<"ERROR, y-bin not found"<<std::endl;
        }
        x_index++;
    }
    std::cout<<"ERROR, x-bin not found"<<std::endl;
}

auto histo_3d::add
(
    std::vector<std::tuple<double, double, double> > news
) noexcept -> void
{
    //Remove all overflow and underflow
    auto uoflow = [&](std::tuple<double, double, double> n)
        {
            if (std::get<0>(n) <= this->xs[0])
            {
                std::get<0>(this->underf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<1>(n) <= this->ys[0])
            {
                std::get<1>(this->underf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<2>(n) <= this->zs[0])
            {
                std::get<2>(this->underf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<0>(n) > this->xs.back())
            {
                std::get<0>(this->overf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<1>(n) > this->ys.back())
            {
                std::get<1>(this->overf)++;
                this->total_counts++;
                return true;
            }
            else if (std::get<2>(n) > this->zs.back())
            {
                std::get<2>(this->overf)++;
                this->total_counts++;
                    return true;
            }
            else
            {
                return false;
            }
        }; 
    std::erase_if(news, uoflow);
    //news.erase
    //(  
    //    std::remove_if(news.begin(), news.end(), uoflow), 
    //    news.end()
    //);

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

    std::vector<std::tuple<double, double> > to_add_xbin{};
    std::vector<double> to_add_xybin{};
    to_add_xbin.reserve(news.size());
    to_add_xybin.reserve(news.size());

    //Sort news by x coordinate
    std::sort
    (
        news.begin(), 
        news.end(),
        [](std::tuple<double, double, double> a, std::tuple<double, double, double> b)
            {
                auto [ a0, ph11, ph12 ] = a;
                auto [ b0, ph21, ph22 ] = b;
                return a0 < b0;
            }
    );
    
    uint16_t x_index = 0, y_index = 0, z_index = 0;
    uint64_t news_index = 0, x_bin_tot = 0, xy_bin_tot = 0,
             to_add_xindex = 0, to_add_xyindex = 0;
    
    //Go through all the x-bins
    for 
    (
        auto x_it = ++this->xs.begin();
        x_it!=this->xs.end(); 
        x_it++
    )
    {
        //Add all the coords that belong in one x-bin at the same time
        while (std::get<0>(news[news_index]) <= *x_it)
        {
            //Stage all coordinates that are in this x-bin
            to_add_xbin.push_back
            (
                std::make_tuple
                (
                    std::get<1>(news[news_index]),
                    std::get<2>(news[news_index])
                )
            );
            news_index++;
            if (news_index == news_tot)
            {
                break;
            }
        }
        if (to_add_xbin.size()>0)
        {
            to_add_xindex = 0;
            x_bin_tot = to_add_xbin.size();

            //Sort the coords by y coordinate
            std::sort
            (
                to_add_xbin.begin(), 
                to_add_xbin.end(),
                [](std::tuple<double, double> a, std::tuple<double, double> b)
                    {
                        auto [ a0, ph1 ] = a;
                        auto [ b0, ph2 ] = b;
                        return a0 < b0;
                    }
            );
            y_index = 0;

            //Go through all the y-bins
            for 
            (
                auto y_it = ++this->ys.begin();
                y_it!=this->ys.end(); 
                y_it++
            )
            {
                //Add all the coords that belong in one y-bin at the same time
                while (std::get<0>(to_add_xbin[to_add_xindex]) <= *y_it)
                {
                    //Stage all coordinates that are in this y-bin
                    to_add_xybin.push_back(std::get<1>(to_add_xbin[to_add_xindex]));
                    to_add_xindex++;
                    if (to_add_xindex == x_bin_tot)
                    {
                        break;
                    }
                }
                if (to_add_xybin.size()>0)
                {
                    to_add_xyindex = 0;
                    xy_bin_tot = to_add_xybin.size();

                    //Sort the coords by z coordinate
                    std::sort(to_add_xybin.begin(), to_add_xybin.end());
                    z_index = 0;

                    //Go through all the z-bins
                    for 
                    (
                        auto z_it = ++this->zs.begin();
                        z_it!=this->zs.end(); 
                        z_it++
                    )
                    {
                        while (to_add_xybin[to_add_xyindex] <= *z_it)
                        {
                            this->counts[x_index][y_index][z_index]++;
                            this->total_counts++;
                            to_add_xyindex++;
                            if (to_add_xyindex == xy_bin_tot)
                            {
                                goto added_staged_xy;
                            }
                        }
                        z_index++;
                    }
                    //All the z-bins checked
                    std::cout<<"ERROR, z-bin not found"<<std::endl;
                }
            added_staged_xy:
                to_add_xybin.clear();
                if (to_add_xindex == x_bin_tot)
                {
                    //All the coords in this x-bin are added
                    goto added_staged;
                }
                y_index++;
            }
            //All the y-bins checked
            std::cout<<"ERROR, y-bin not found"<<std::endl;
        }
    added_staged:
        if (news_index == news_tot)
        {
            //All the new coords are added
            return;
        }
        to_add_xbin.clear();
        x_index++;
    }
    //All the x-bins checked
    std::cout<<"ERROR, x-bin not found"<<std::endl;
}

auto histo_3d::add
(
    const histo_3d &other
) noexcept -> void
{
    if (!(*this == other))
    {
        std::cout<<"Incompatible histograms!"<<std::endl;
        return;
    }

    std::get<0>(this->underf) += std::get<0>(other.underf);
    std::get<1>(this->underf) += std::get<1>(other.underf);
    std::get<2>(this->underf) += std::get<2>(other.underf);

    std::get<0>(this->overf) += std::get<0>(other.overf);
    std::get<1>(this->overf) += std::get<1>(other.overf);
    std::get<2>(this->overf) += std::get<2>(other.overf);

    for (uint16_t i = 0; i<this->counts.size(); i++)
    {
        for (uint16_t j = 0; j<this->counts[0].size(); j++)
        {
            for (uint16_t k = 0; k<this->counts[0][0].size(); k++)
            {
                this->counts[i][j][k] += other.counts[i][j][k];
            }
        }
    }

    this->total_counts += other.total_counts;
}

auto histo_1d::get_histo() const noexcept 
-> std::tuple<std::vector<double>, double, double, std::vector<double>, uint64_t>
{   
    auto ret_histo{this->counts};
    auto total = this->total_counts;

    auto prev_x = this->xs[0];
    double curr_x, bin_size;

    for (uint16_t i = 0; i<this->counts.size(); i++)
    {
        curr_x = this->xs[i+1];
        bin_size = curr_x-prev_x;
        prev_x = curr_x;
        ret_histo[i] /= static_cast<double>(total)*bin_size;

    }

    return std::make_tuple
        (
            ret_histo,
            this->underf / static_cast<double>(total),
            this->overf / static_cast<double>(total),
            this->xs, 
            this->total_counts
        );
}

auto histo_2d::get_histo() const noexcept -> std::tuple
<   std::vector<std::vector<double> >, 
    std::tuple<double, double>, 
    std::tuple<double, double>, 
    std::tuple<std::vector<double>, std::vector<double> >,
    uint64_t >
{   
    auto ret_histo{this->counts};
    auto total = this->total_counts;

    auto prev_x = this->xs[0];
    double curr_x, x_bin_size;
    double prev_y, curr_y, y_bin_size;

    for (uint16_t i = 0; i<this->counts.size(); i++)
    {
        curr_x = this->xs[i+1];
        x_bin_size = curr_x-prev_x;
        prev_x = curr_x;

        prev_y = this->ys[0];

        for (uint16_t j = 0; j<this->counts[0].size(); j++)
        {
            curr_y = this->ys[j+1];
            y_bin_size = curr_y-prev_y;
            prev_y = curr_y;

            ret_histo[i][j] /= static_cast<double>(total)*x_bin_size*y_bin_size;
        }
    }

    return std::make_tuple
        (
            std::move(ret_histo),
            std::make_tuple(std::get<0>(this->underf) / static_cast<double>(total), 
                            std::get<1>(this->underf) / static_cast<double>(total)),
            std::make_tuple(std::get<0>(this->overf) / static_cast<double>(total), 
                            std::get<1>(this->overf) / static_cast<double>(total)),
            std::make_tuple(this->xs, 
                            this->ys),
            this->total_counts
        );
}

auto histo_3d::get_histo() const noexcept -> std::tuple
<   std::vector<std::vector<std::vector<double> > >, 
    std::tuple<double, double, double>, 
    std::tuple<double, double, double>, 
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double> >,
    uint64_t >
{   
    auto ret_histo{this->counts};
    auto total = this->total_counts;

    auto n_x_bins = this->counts.size();
    auto n_y_bins = this->counts[0].size();
    auto n_z_bins = this->counts[0][0].size();
    
    std::vector<double> x_bins(n_x_bins);
    std::vector<double> y_bins(n_y_bins);
    std::vector<double> z_bins(n_z_bins);

    {
        double prev, curr;
        prev = this->xs[0];
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            curr = this->xs[i+1];
            x_bins[i] = curr - prev;
            prev = curr;
        }
        prev = this->ys[0];
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            curr = this->ys[i+1];
            y_bins[i] = curr - prev;
            prev = curr;
        }
        prev = this->zs[0];
        for (uint16_t i = 0; i < n_z_bins; i++)
        {
            curr = this->zs[i+1];
            z_bins[i] = curr - prev;
            prev = curr;
        }
    }

    for (uint16_t i = 0; i < n_x_bins; i++)
    {
        for (uint16_t j = 0; j < n_y_bins; j++)
        {
            for (uint16_t k = 0; k < n_z_bins; k++)
            {
                ret_histo[i][j][k] /= static_cast<double>(total)*x_bins[i]*y_bins[j]*z_bins[k];
            }
        }
    }

    return std::make_tuple
        (
            std::move(ret_histo),
            std::make_tuple(std::get<0>(this->underf) / static_cast<double>(total), 
                            std::get<1>(this->underf) / static_cast<double>(total), 
                            std::get<2>(this->underf) / static_cast<double>(total)),
            std::make_tuple(std::get<0>(this->overf) / static_cast<double>(total), 
                            std::get<1>(this->overf) / static_cast<double>(total), 
                            std::get<2>(this->overf) / static_cast<double>(total)),
            std::make_tuple(this->xs, 
                            this->ys, 
                            this->zs),
            this->total_counts
        );
}

auto histo_2d::project_1d(const bool project_ys) const noexcept -> std::tuple
<   std::vector<double>&&, 
    double, 
    double, 
    std::vector<double>&&, 
    uint64_t >
{   
    auto [ full_histo, uf, of, xs, tots ] = this->get_histo();
    std::vector<double> ret_histo;
    std::vector<double> ret_grid;
    std::vector<double> bins;
    double curr, prev;
    uint16_t n_x_bins = static_cast<uint16_t>(std::get<0>(xs).size()-1);
    uint16_t n_y_bins = static_cast<uint16_t>(std::get<1>(xs).size()-1);

    if (project_ys == true)
    {
        bins.reserve(n_y_bins);
        prev = std::get<1>(xs)[0];
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            curr = std::get<1>(xs)[i+1];
            bins[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<double>(n_x_bins, 0.0);
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            for (uint16_t j = 0; j < n_y_bins; j++)
            {
                ret_histo[i] += full_histo[i][j]*bins[j];
            }
        }

        ret_grid = std::get<0>(xs);
    }
    else
    {
        bins.reserve(n_x_bins);
        prev = std::get<0>(xs)[0];
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            curr = std::get<0>(xs)[i+1];
            bins[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<double>(n_y_bins, 0.0);
        for (uint16_t i = 0; i<n_y_bins; i++)
        {
            for (uint16_t j = 0; j<n_x_bins; j++)
            {
                ret_histo[i] += full_histo[j][i]*bins[j];
            }
        }

        ret_grid = std::get<1>(xs);
    }

    return std::make_tuple
        (
            std::move(ret_histo),
            std::get<0>(uf)+std::get<1>(uf),
            std::get<0>(of)+std::get<1>(of),
            std::move(ret_grid), 
            tots
        );
}

auto histo_3d::project_1d(const uint8_t dim_left) const noexcept -> std::tuple
<   std::vector<double>&&, 
    double, 
    double, 
    std::vector<double>&&, 
    uint64_t >
{   
    auto [ full_histo, uf, of, xs, tots ] = this->get_histo();
    std::vector<double> ret_histo;
    std::vector<double> ret_grid;
    std::vector<double> bins1;
    std::vector<double> bins2;
    double curr, prev;
    uint16_t n_x_bins = static_cast<uint16_t>(std::get<0>(xs).size()-1);
    uint16_t n_y_bins = static_cast<uint16_t>(std::get<1>(xs).size()-1);
    uint16_t n_z_bins = static_cast<uint16_t>(std::get<2>(xs).size()-1);

    if (dim_left == 0) //x
    {
        bins1.reserve(n_y_bins); //y
        bins2.reserve(n_z_bins); //z

        prev = std::get<1>(xs)[0];
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            curr = std::get<1>(xs)[i+1];
            bins1[i] = curr - prev;
            prev = curr;
        }

        prev = std::get<2>(xs)[0];
        for (uint16_t i = 0; i < n_z_bins; i++)
        {
            curr = std::get<2>(xs)[i+1];
            bins2[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<double>(n_x_bins, 0.0);
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            for (uint16_t j = 0; j < n_y_bins; j++)
            {
                for (uint16_t k = 0; k < n_z_bins; k++)
                {
                    ret_histo[i] += full_histo[i][j][k]*bins1[j]*bins2[k];
                }
            }
        }

        ret_grid = std::get<0>(xs);
    }
    else if (dim_left == 1) //y
    {
        bins1.reserve(n_x_bins); //x
        bins2.reserve(n_z_bins); //z

        prev = std::get<0>(xs)[0];
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            curr = std::get<0>(xs)[i+1];
            bins1[i] = curr - prev;
            prev = curr;
        }
        
        prev = std::get<2>(xs)[0];
        for (uint16_t i = 0; i < n_z_bins; i++)
        {
            curr = std::get<2>(xs)[i+1];
            bins2[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<double>(n_y_bins, 0.0);
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            for (uint16_t j = 0; j < n_x_bins; j++)
            {
                for (uint16_t k = 0; k < n_z_bins; k++)
                {
                    ret_histo[i] += full_histo[j][i][k]*bins1[j]*bins2[k];
                }
            }
        }

        ret_grid = std::get<1>(xs);
    }
    else if (dim_left == 2) //z
    {
        bins1.reserve(n_y_bins); //y
        bins2.reserve(n_z_bins); //x

        prev = std::get<1>(xs)[0];
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            curr = std::get<1>(xs)[i+1];
            bins1[i] = curr - prev;
            prev = curr;
        }
        
        prev = std::get<0>(xs)[0];
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            curr = std::get<0>(xs)[i+1];
            bins2[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<double>(n_z_bins, 0.0);
        for (uint16_t i = 0; i < n_z_bins; i++)
        {
            for (uint16_t j = 0; j < n_x_bins; j++)
            {
                for (uint16_t k = 0; k < n_y_bins; k++)
                {
                    ret_histo[i] += full_histo[j][k][i]*bins1[k]*bins2[j];
                }
            }
        }

        ret_grid = std::get<2>(xs);
    }
    else
    {
        std::cout<<"ERROR dimension must be 0=x, 1=y or 2=z, given: ";
        std::cout<<dim_left<<std::endl;
        std::cout<<"Output will be garbage"<<std::endl;
    }

    return std::make_tuple
        (
            std::move(ret_histo),
            std::get<0>(uf)+std::get<1>(uf)+std::get<2>(uf),
            std::get<0>(of)+std::get<1>(of)+std::get<2>(of),
            std::move(ret_grid), 
            tots
        );
}
auto histo_3d::project_2d(const uint8_t dim_integrated)
const noexcept -> std::tuple
<   std::vector<std::vector<double> >, 
    double, 
    double, 
    std::tuple<std::vector<double>, std::vector<double> >,
    uint64_t >
{
    auto [ full_histo, uf, of, xs, tots ] = this->get_histo();
    std::vector<std::vector<double> > ret_histo;
    std::vector<double> ret_grid1;
    std::vector<double> ret_grid2;
    std::vector<double> bins;
    double curr, prev;
    uint16_t n_x_bins = static_cast<uint16_t>(std::get<0>(xs).size()-1);
    uint16_t n_y_bins = static_cast<uint16_t>(std::get<1>(xs).size()-1);
    uint16_t n_z_bins = static_cast<uint16_t>(std::get<2>(xs).size()-1);

    if (dim_integrated == 0) //x
    {
        bins = std::vector<double>(n_x_bins, 0.0);

        prev = std::get<0>(xs)[0];
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            curr = std::get<0>(xs)[i+1];
            bins[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<std::vector<double> >
                    (
                        n_y_bins, 
                        std::vector<double>
                        (
                            n_z_bins, 
                            0.0
                        )
                    );
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            for (uint16_t j = 0; j < n_z_bins; j++)
            {
                for (uint16_t k = 0; k < n_x_bins; k++)
                {
                    ret_histo[i][j] += full_histo[k][i][j]*bins[k];
                }
            }
        }

        ret_grid1 = std::get<1>(xs);
        ret_grid2 = std::get<2>(xs);
    }
    else if (dim_integrated == 1) //y
    {
        bins = std::vector<double>(n_y_bins, 0.0);

        prev = std::get<1>(xs)[0];
        for (uint16_t i = 0; i < n_y_bins; i++)
        {
            curr = std::get<1>(xs)[i+1];
            bins[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<std::vector<double> >
                    (
                        n_x_bins, 
                        std::vector<double>
                        (
                            n_z_bins, 
                            0.0
                        )
                    );
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            for (uint16_t j = 0; j < n_z_bins; j++)
            {
                for (uint16_t k = 0; k < n_y_bins; k++)
                {
                    ret_histo[i][j] += full_histo[i][k][j]*bins[k];
                }
            }
        }

        ret_grid1 = std::get<0>(xs);
        ret_grid2 = std::get<2>(xs);
    }
    else if (dim_integrated == 2) //z
    {
        bins = std::vector<double>(n_z_bins, 0.0);

        prev = std::get<2>(xs)[0];
        for (uint16_t i = 0; i < n_z_bins; i++)
        {
            curr = std::get<2>(xs)[i+1];
            bins[i] = curr - prev;
            prev = curr;
        }

        ret_histo = std::vector<std::vector<double> >
                    (
                        n_x_bins, 
                        std::vector<double>
                        (
                            n_y_bins, 
                            0.0
                        )
                    );
        for (uint16_t i = 0; i < n_x_bins; i++)
        {
            for (uint16_t j = 0; j < n_y_bins; j++)
            {
                for (uint16_t k = 0; k < n_z_bins; k++)
                {
                    ret_histo[i][j] += full_histo[i][j][k]*bins[k];
                }
            }
        }

        ret_grid1 = std::get<0>(xs);
        ret_grid2 = std::get<1>(xs);
    }
    else
    {
        std::cout<<"ERROR dimension must be 0=x, 1=y or 2=z, given: ";
        std::cout<<dim_integrated<<std::endl;
        std::cout<<"Output will be garbage"<<std::endl;
    }

    return std::make_tuple
        (
            std::move(ret_histo),
            std::get<0>(uf)+std::get<1>(uf)+std::get<2>(uf),
            std::get<0>(of)+std::get<1>(of)+std::get<2>(of),
            std::make_tuple(ret_grid1, ret_grid2), 
            tots
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

auto operator==
(
    const histo_2d& c1, 
    const histo_2d& c2
) noexcept -> bool
{
    bool ret = (c1.xs.size() == c2.xs.size());
    ret *= (c1.ys.size() == c2.ys.size());
    for (uint16_t i = 0; i<c1.xs.size(); i++)
    {
        ret *= (c1.xs[i] == c2.xs[i]);
    }
    for (uint16_t i = 0; i<c1.ys.size(); i++)
    {
        ret *= (c1.ys[i] == c2.ys[i]);
    }
    return ret;
}

auto operator==
(
    const histo_3d& c1, 
    const histo_3d& c2
) noexcept -> bool
{
    bool ret = (c1.xs.size() == c2.xs.size());
    ret *= (c1.ys.size() == c2.ys.size());
    ret *= (c1.zs.size() == c2.zs.size());
    for (uint16_t i = 0; i<c1.xs.size(); i++)
    {
        ret *= (c1.xs[i] == c2.xs[i]);
    }
    for (uint16_t i = 0; i<c1.ys.size(); i++)
    {
        ret *= (c1.ys[i] == c2.ys[i]);
    }
    for (uint16_t i = 0; i<c1.zs.size(); i++)
    {
        ret *= (c1.zs[i] == c2.zs[i]);
    }
    return ret;
}
/*
int main()
{
    std::vector<double> xs{0.0,0.5,1.0,1.5,2.5};
    histo_1d h1{xs}, h2{xs};
    auto [ c1, uf1, of1, x1, tot1 ] = h1.get_histo();
    
    for (auto x : c1) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto x : x1) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf1<<std::endl;
    std::cout<<of1<<std::endl;
    std::cout<<tot1<<std::endl;
    std::cout<<std::endl;

    h1.add(0.1);
    h1.add(std::vector<double>({1.1,0.2,-2}));
    auto [ c2, uf2, of2, x2, tot2 ] = h1.get_histo();
    
    for (auto x : c2) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto x : x2) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf2<<std::endl;
    std::cout<<of2<<std::endl;
    std::cout<<tot2<<std::endl;
    std::cout<<std::endl;

    h2.add(2.0);
    h1.add(h2);
    auto [ c3, uf3, of3, x3, tot3 ] = h1.get_histo();
    
    for (auto x : c3) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto x : x3) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf3<<std::endl;
    std::cout<<of3<<std::endl;
    std::cout<<tot3<<std::endl;
    std::cout<<std::endl;

    std::vector<double> ys{1.0,1.5,2.5,3.0};
    histo_2d h2d{xs, ys};
    auto [ c4, uf4, of4, x4, tot4 ] = h2d.get_histo();
    auto [ uf40, uf41 ] = uf4;
    auto [ of40, of41 ] = of4;
    auto [ x40, x41 ] = x4;
    
    for (auto x : c4)
    {
        for (auto y : x)
            std::cout<<y<<' ';
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    for (auto x : x40) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto y : x41) std::cout<<y<<' ';
    std::cout<<std::endl;
    std::cout<<uf40<<' '<<uf41<<std::endl;
    std::cout<<of40<<' '<<of41<<std::endl;
    std::cout<<tot4<<std::endl;
    std::cout<<std::endl;

    h2d.add({0.1, 0.2});
    h2d.add(std::vector<std::tuple<double, double> >({{1.1,0.2},{1.1,1.2},{1.1,2.2},{2.1,2.7},{4.1,0.2}}));
    auto [ c5, uf5, of5, x5, tot5 ] = h2d.get_histo();
    auto [ uf50, uf51 ] = uf5;
    auto [ of50, of51 ] = of5;
    auto [ x50, x51 ] = x5;
    
    for (auto x : c5)
    {
        for (auto y : x)
            std::cout<<y<<' ';
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    for (auto x : x50) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto y : x51) std::cout<<y<<' ';
    std::cout<<std::endl;
    std::cout<<uf50<<' '<<uf51<<std::endl;
    std::cout<<of50<<' '<<of51<<std::endl;
    std::cout<<tot5<<std::endl;
    std::cout<<std::endl;

    histo_2d h2d2{xs, ys};
    h2d2.add({-0.1, 0.2});
    h2d2.add(std::vector<std::tuple<double, double> >({{1.1,1.2},{1.2,2.2},{1.1,2.5},{0.1,2.7},{0.1,1.2}}));
    h2d2.add(std::vector<std::tuple<double, double> >({{0.1,0.2},{4.1,1.2},{1.1,5.2},{0.1,1.7}}));
    h2d.add(h2d2);
    auto [ c6, uf6, of6, x6, tot6 ] = h2d.get_histo();
    auto [ uf60, uf61 ] = uf6;
    auto [ of60, of61 ] = of6;
    auto [ x60, x61 ] = x6;
    
    for (auto x : c6)
    {
        for (auto y : x)
            std::cout<<y<<' ';
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    for (auto x : x60) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto y : x61) std::cout<<y<<' ';
    std::cout<<std::endl;
    std::cout<<uf60<<' '<<uf61<<std::endl;
    std::cout<<of60<<' '<<of61<<std::endl;
    std::cout<<tot6<<std::endl;
    std::cout<<std::endl;
    
    auto [ c7, uf7, of7, x7, tot7 ] = h2d.project_1d(true);
    
    for (auto x : c7)
    {
        std::cout<<x<<' ';
    }
    std::cout<<std::endl;
    for (auto x : x7) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf7<<std::endl;
    std::cout<<of7<<std::endl;
    std::cout<<tot7<<std::endl;
    std::cout<<std::endl;
    
    auto [ c8, uf8, of8, x8, tot8 ] = h2d.project_1d(false);
    
    for (auto x : c8)
    {
        std::cout<<x<<' ';
    }
    std::cout<<std::endl;
    for (auto x : x8) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf8<<std::endl;
    std::cout<<of8<<std::endl;
    std::cout<<tot8<<std::endl;
    std::cout<<std::endl;

    xs = std::vector<double>{0, 1, 1.5};
    ys = std::vector<double>{1, 2, 2.5};
    auto zs = std::vector<double>{2, 3, 3.5};

    histo_3d h3d{xs, ys, zs};
    h3d.add({0.5, 2.2, 3.2});
    h3d.add(std::vector<std::tuple<double, double, double> >
        ({{-1.1, 1.2, 3},{2, 2.2, 3},{1.1, 0.5, 3},{1.1, 2.7, 3},{1.1, 1.2, 1},{1.1, 1.2, 4}}));
    h3d.add(std::vector<std::tuple<double, double, double> >
        ({{1.1, 1.2, 3},{2, 2.2, 2.3},{0.1, 2.3, 3.1},{1.1, 1.7, 2.3},{0.2, 1.3, 3.1},{1.1, 1.2, 3.3}}));
    
    auto [ c9, uf9, of9, x9, tot9 ] = h3d.get_histo();
    auto [ uf90, uf91, uf92 ] = uf9;
    auto [ of90, of91, of92 ] = of9;
    auto [ x90, x91, x92 ] = x9;
    
    std::cout<<std::endl;
    for (auto x : c9)
    {
        for (auto y : x)
        {
            for (auto z : y)
                std::cout<<z<<' ';
            std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    for (auto x : x90) std::cout<<x<<' ';
    std::cout<<std::endl;
    for (auto y : x91) std::cout<<y<<' ';
    std::cout<<std::endl;
    for (auto z : x92) std::cout<<z<<' ';
    std::cout<<std::endl;
    std::cout<<uf90<<' '<<uf91<<' '<<uf92<<std::endl;
    std::cout<<of90<<' '<<of91<<' '<<of92<<std::endl;
    std::cout<<tot9<<std::endl;
    std::cout<<std::endl;
    
    auto [ c10, uf10, of10, x10, tot10 ] = h3d.project_1d(0);
    
    for (auto x : c10)
    {
        std::cout<<x<<' ';
    }
    std::cout<<std::endl;
    for (auto x : x10) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf10<<std::endl;
    std::cout<<of10<<std::endl;
    std::cout<<tot10<<std::endl;
    std::cout<<std::endl;
    
    auto [ c11, uf11, of11, x11, tot11 ] = h3d.project_1d(1);
    
    for (auto x : c11)
    {
        std::cout<<x<<' ';
    }
    std::cout<<std::endl;
    for (auto x : x11) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf11<<std::endl;
    std::cout<<of11<<std::endl;
    std::cout<<tot11<<std::endl;
    std::cout<<std::endl;
    
    auto [ c12, uf12, of12, x12, tot12 ] = h3d.project_1d(2);
    
    for (auto x : c12)
    {
        std::cout<<x<<' ';
    }
    std::cout<<std::endl;
    for (auto x : x12) std::cout<<x<<' ';
    std::cout<<std::endl;
    std::cout<<uf12<<std::endl;
    std::cout<<of12<<std::endl;
    std::cout<<tot12<<std::endl;
    std::cout<<std::endl;

    return 0;
}*/