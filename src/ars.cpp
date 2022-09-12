//Copyright (c) 2021 Mikko Kuha

#include "ars.hpp"

ars::ars(const std::function<double(const double&)> & pdf_, const double & min, const double & max) noexcept 
    : adaptive{true}, absc_count_limit{20}, log_pdf([&](const double& x){return log(pdf_(x));}), log_der([&](const double& x){return (log(pdf_(x+1e-5))-log(pdf_(x-1e-5)))/2e-5;}), abscissae({min+(max-min)/10.0, max-(max-min)/10.0}), zk({min,max}), uk(), unif_dist(0.0,1.0)
{
    std::unique_lock lock{this->m};
    this->log_values.clear();
    this->log_der_values.clear();

    this->initialize_distribution();
}

ars::ars(const std::function<double(const double&)> & pdf_, const std::function<double(const double&)> & pdf_derivative_, const std::initializer_list<double> & init_abscissae, const std::initializer_list<double> & init_intersections, const std::initializer_list<double> & init_uk, const std::map<double, double> & init_log_values, const std::map<double, double> & init_log_der_values) noexcept
    : adaptive{true}, absc_count_limit{20}, log_pdf([&](const double& x){return log(pdf_(x));}), log_der([&](const double& x){return log(pdf_derivative_(x));}), log_values(init_log_values), log_der_values(init_log_der_values), abscissae(init_abscissae), zk(init_intersections), uk(init_uk), unif_dist(0.0,1.0) 
{
    std::unique_lock lock{this->m};
    auto temp = std::vector<double>(this->uk.size(), 0.0);
    auto uk_iter = this->uk.cbegin();
    auto temp_iter = temp.begin();
    auto temp_end = temp.cend();
    for (; temp_iter!=temp_end; temp_iter++, uk_iter++)
    {
        *temp_iter = exp(*uk_iter);
    }

    //this->report();

    this->sk = ars::piecewise_exponential_distribution(this->zk.cbegin(), this->zk.cend(), temp.cbegin());
}

double ars::throw_one_adaptive(std::mt19937 & random_generator)
{
    std::unique_lock lock{this->m};

    if (!(this->abscissae.size()<this->absc_count_limit))
    {
        this->adaptive = false;
        lock.~unique_lock();
        return this->throw_one_const(random_generator);
    }
    do
    {
        try
        {
            const double rand = this->sk(random_generator);
            const double comparator = this->unif_dist(random_generator);

            //Index of the lower bound of the zk interval where rand falls
            uint k=0; 
            for (; k<this->zk.size()-1; k++)
            {
                if (this->zk.at(k+1) > rand) break;
            }
            const double xk = this->abscissae.at(k);
            //Index of the lower bound of the xk interval where rand falls
            const int temp_j = (xk > rand)? static_cast<int>(k-1) : static_cast<int>(k);
            //rand < x0  OR > xlast -> re-sample
            if (temp_j==-1 || temp_j==static_cast<int>(this->abscissae.size())-1) 
            {
                if (this->abscissae.size()<this->absc_count_limit)
                {
                    this->abscissae.push_back(rand);
                    std::sort(this->abscissae.begin(),this->abscissae.end());
                    this->update_distribution();
                }//std::cout<<"HERE 1"<<rand<<std::endl;
                continue;
            }
            const uint j = static_cast<uint>(temp_j);

            const double hxk = this->log_values.at(xk);
            const double hdxk = this->log_der_values.at(xk);
            const double ukx = hxk + (rand - xk)*hdxk;

            double ljx=0, upper_hxj=0, lower_hxj=0, upper_xj=0, lower_xj=0;
            if (j==k)
            {
                lower_xj = xk; 
                upper_xj = this->abscissae.at(k+1);
            }
            else
            {
                lower_xj = this->abscissae.at(k-1); 
                upper_xj = xk;
            }
            lower_hxj = this->log_values.at(lower_xj);
            upper_hxj = this->log_values.at(upper_xj);

            ljx = ((upper_xj-rand)*lower_hxj + (rand-lower_xj)*upper_hxj)/(upper_xj - lower_xj);

            if (comparator <= exp(ljx - ukx))
            {//std::cout<<"HERE 2 "<<rand<<std::endl;
                return rand;
            }
            if (this->abscissae.size()<this->absc_count_limit)
            {
                this->abscissae.push_back(rand);
                std::sort(this->abscissae.begin(),this->abscissae.end());
                this->update_distribution();
            
                if (comparator <= exp(this->log_values.at(rand) - ukx))
                {//std::cout<<"HERE 3 "<<rand<<std::endl;
                    return rand;
                }
            }
            else if (comparator <= exp(this->log_pdf(rand) - ukx))
            {//std::cout<<"HERE 4 "<<rand<<std::endl;
                return rand;
            }

        }
        catch(const std::exception& e)
        {
            std::cout << e.what() << std::endl << std::endl;
        }
    } while (true);
}   

double ars::throw_one_const(std::mt19937 & random_generator)
{
    do
    {
        try
        {
            const double rand = this->sk(random_generator);
            const double comparator = this->unif_dist(random_generator);

            //Index of the lower bound of the zk interval where rand falls
            uint k=0; 
            for (; k<this->zk.size()-1; k++)
            {
                if (this->zk.at(k+1) > rand) break;
            }
            const double xk = this->abscissae.at(k);
            //Index of the lower bound of the xk interval where rand falls
            const int temp_j = (xk > rand)? static_cast<int>(k-1) : static_cast<int>(k);
            //rand < x0  OR > xlast -> re-sample
            if (temp_j==-1 || temp_j==static_cast<int>(this->abscissae.size())-1) 
            {//std::cout<<"HERE 1"<<rand<<std::endl;
                continue;
            }
            const uint j = static_cast<uint>(temp_j);

            const double hxk = this->log_values.at(xk);
            const double hdxk = this->log_der_values.at(xk);
            const double ukx = hxk + (rand - xk)*hdxk;

            double ljx=0, upper_hxj=0, lower_hxj=0, upper_xj=0, lower_xj=0;
            if (j==k)
            {
                lower_xj = xk; 
                upper_xj = this->abscissae.at(k+1);
            }
            else
            {
                lower_xj = this->abscissae.at(k-1); 
                upper_xj = xk;
            }
            lower_hxj = this->log_values.at(lower_xj);
            upper_hxj = this->log_values.at(upper_xj);

            ljx = ((upper_xj-rand)*lower_hxj + (rand-lower_xj)*upper_hxj)/(upper_xj - lower_xj);

            if (comparator <= exp(ljx - ukx))
            {//std::cout<<"HERE 2 "<<rand<<std::endl;
                return rand;
            }
            if (comparator <= exp(this->log_pdf(rand) - ukx))
            {//std::cout<<"HERE 4 "<<rand<<std::endl;
                return rand;
            }
        }
        catch(const std::exception& e)
        {
            std::cout << e.what() << " in ars, trying again"<<std::endl;
        }
    } while (true);
}   

std::vector<double> ars::throw_n(const uint & n, std::mt19937 & random_generator) noexcept
{
    std::vector<double> vec({});
    vec.reserve(n); 
    for (uint i=0; i<n; i++)
    {
        vec.emplace_back(this->throw_one(random_generator));
    }
    return vec;
}


void ars::initialize_distribution() noexcept
{
    for (uint k=0; k<this->abscissae.size(); k++)
    {
        auto xk = this->abscissae.at(k);

        this->log_values.emplace(xk,this->log_pdf(xk));
        this->log_der_values.emplace(xk,this->log_der(xk));
    }

    for (uint k=0; k<this->abscissae.size()-1; k++)
    {
        auto xk = this->abscissae.at(k);
        auto xkp1 = this->abscissae.at(k+1);
        auto iter = this->log_values.find(xk);
        auto iter2 = this->log_der_values.find(xk);

        double hk, hkp1, hdk, hdkp1;
        std::tie(std::ignore, hk) = *iter;
        std::tie(std::ignore, hkp1) = *std::next(iter);
        std::tie(std::ignore, hdk) = *iter2;
        std::tie(std::ignore, hdkp1) = *std::next(iter2);

        this->zk.emplace(std::next(this->zk.cbegin()), (hkp1 - hk - xkp1*hdkp1 + xk*hdk)/(hdk - hdkp1));
    }

    this->uk = std::vector<double>(zk.size(), 0.0);
    auto weights_iter = this->uk.begin();
    auto zj_iter = this->zk.cbegin();
    auto hxj_iter = this->log_values.cbegin();
    auto hdxj_iter = this->log_der_values.cbegin();
    auto endvalue = this->log_values.cend();
    for (; hxj_iter!=endvalue; weights_iter++, zj_iter++, hxj_iter++, hdxj_iter++)
    {
        *weights_iter = (*hxj_iter).second + (*zj_iter - (*hxj_iter).first)*((*hdxj_iter).second);
    }
    this->uk.back() = (*this->log_values.crbegin()).second + (*this->zk.crbegin() - (*this->log_values.crbegin()).first)*((*this->log_der_values.crbegin()).second);
    
    auto temp = std::vector<double>(this->uk.size(), 0.0);
    auto uk_iter = this->uk.cbegin();
    auto temp_iter = temp.begin();
    auto temp_end = temp.cend();
    for (; temp_iter!=temp_end; temp_iter++, uk_iter++)
    {
        *temp_iter = exp(*uk_iter);
    }

    //this->report();

    this->sk = ars::piecewise_exponential_distribution(this->zk.cbegin(), this->zk.cend(), temp.cbegin());
}

void ars::update_distribution() noexcept
{
    bool new_abscissae = false;

    //std::cout<<"update\n";

    auto x0 = this->abscissae.front();
    auto iter0 = this->log_values.find(x0);
    if (iter0 == this->log_values.cend())
    {
        new_abscissae = true;
        double hk, hkp1, hdk, hdkp1;
        auto xkp1 = this->abscissae.at(1);
        auto temp = this->log_values.emplace(x0,this->log_pdf(x0));
        std::tie(std::ignore, hk) = *temp.first;
        std::tie(std::ignore, hkp1) = *std::next(temp.first);

        auto temp2 = this->log_der_values.emplace(x0,this->log_der(x0));
        std::tie(std::ignore, hdk) = *temp2.first;
        std::tie(std::ignore, hdkp1) = *std::next(temp2.first);
        this->zk.emplace(std::next(this->zk.cbegin()), (hkp1 - hk - xkp1*hdkp1 + x0*hdk)/(hdk - hdkp1));
    }

    for (uint k=1; k<this->abscissae.size()-1; k++)
    {
        auto xk = this->abscissae.at(k);
        auto iter = this->log_values.find(xk);
        if (iter == this->log_values.cend())
        {
            new_abscissae = true;
            double hk, hkm1, hdk, hdkm1;
            auto xkm1 = this->abscissae.at(k-1);

            auto temp = this->log_values.emplace(xk,this->log_pdf(xk));
            std::tie(std::ignore, hk) = *temp.first;
            std::tie(std::ignore, hkm1) = *std::prev(temp.first);
            
            auto temp2 = this->log_der_values.emplace(xk,this->log_der(xk));
            std::tie(std::ignore, hdk) = *temp2.first;
            std::tie(std::ignore, hdkm1) = *std::prev(temp2.first);

            auto prev_value = this->zk.at(k);
            this->zk.at(k) = (hk - hkm1 - xk*hdk + xkm1*hdkm1)/(hdkm1 - hdk);

            if (k != this->abscissae.size()-1)
            {
                double hkp1, hdkp1;
                auto xkp1 = this->abscissae.at(k+1);
                std::tie(std::ignore, hkp1) = *std::next(temp.first);
                std::tie(std::ignore, hdkp1) = *std::next(temp2.first);
                this->zk.emplace(std::next(this->zk.cbegin(),k+1), (hkp1 - hk - xkp1*hdkp1 + xk*hdk)/(hdk - hdkp1));
            }
            else
            {
                this->zk.emplace_back(prev_value);
            }            
        }
    }

    if (new_abscissae)
    {
        this->uk = std::vector<double>(this->zk.size(), 0.0);
        auto weights_iter = this->uk.begin();
        auto zj_iter = this->zk.cbegin();
        auto hxj_iter = this->log_values.cbegin();
        auto hdxj_iter = this->log_der_values.cbegin();
        auto endvalue = this->log_values.cend();
        for (; hxj_iter!=endvalue; weights_iter++, zj_iter++, hxj_iter++, hdxj_iter++)
        {
            *weights_iter = (*hxj_iter).second + (*zj_iter - (*hxj_iter).first)*((*hdxj_iter).second);
        }
        this->uk.back() = (*this->log_values.crbegin()).second + (*this->zk.crbegin() - (*this->log_values.crbegin()).first)*((*this->log_der_values.crbegin()).second);

        auto temp = std::vector<double>(this->uk.size(), 0.0);
        auto uk_iter = this->uk.cbegin();
        auto temp_iter = temp.begin();
        auto temp_end = temp.cend();
        for (; temp_iter!=temp_end; temp_iter++, uk_iter++)
        {
            *temp_iter = exp(*uk_iter);
        }
        //this->report();

        this->sk = ars::piecewise_exponential_distribution(this->zk.cbegin(), this->zk.cend(), temp.cbegin());
    }
}