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
        for (auto & r : ret_histo)
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

    auto histo_2d::get_histo() const noexcept
    {   
        std::vector<std::vector<double> > ret_histo{this->counts};
        auto total = this->total_counts;
        for (auto & row : ret_histo)
        {
            for (auto & r : row)
            {
                r /= total;
            }
        }

        return std::make_tuple
            (
                std::move(ret_histo),
                std::make_tuple(std::get<0>(this->underf) / total, 
                                std::get<1>(this->underf) / total),
                std::make_tuple(std::get<0>(this->overf) / total, 
                                std::get<1>(this->overf) / total),
                std::make_tuple(this->xs, 
                                this->ys),
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

        return 0;
    }