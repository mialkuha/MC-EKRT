// Copyright (c) 2025 Mikko Kuha (University of Jyväskylä).
// This program is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details. You should have received a copy
// of the GNU General Public License along with this program. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef ARS_HPP
#define ARS_HPP

#include <iostream>

#include <cmath>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <functional>
#include <random>
#include <vector>

#include "typedefs.hpp"

class ars
{
public:
    class piecewise_exponential_distribution
    {
    public:
        piecewise_exponential_distribution() : normalization{}, abscissae{}, weights{} {};
        piecewise_exponential_distribution &operator=(const piecewise_exponential_distribution &other)
        {
            if (this != &other)
            {
                std::unique_lock lock{this->m};
                this->normalization = other.normalization;
                this->exp_factors = other.exp_factors;
                this->cdf_weights = other.cdf_weights;
                this->abscissae = other.abscissae;
                this->weights = other.weights;
            }
            return *this;
        }
        template <class InputIteratorB, class InputIteratorW>
        piecewise_exponential_distribution(InputIteratorB firstB, InputIteratorB lastB, InputIteratorW firstW) : normalization{}, exp_factors{}, cdf_weights{}, abscissae{}, weights{}
        {
            std::unique_lock lock{this->m};
            for (; firstB != lastB; firstB++, firstW++)
            {
                this->abscissae.emplace_back(*firstB);
                this->weights.emplace_back(*firstW);
            }

            this->initialize_members();
        }
        void add_point(const double &absc, const double &weight) noexcept
        {
            std::unique_lock lock{this->m};
            uint_fast16_t i = 0;
            for (; i < this->abscissae.size(); i++)
            {
                if (absc < this->abscissae.at(i))
                {
                    break;
                }
            }
            this->abscissae.insert(this->abscissae.begin() + static_cast<int_fast32_t>(i), absc);
            this->weights.insert(this->weights.begin() + static_cast<int_fast32_t>(i), weight);

            this->initialize_members();
            /* maybe TODO to be more efficient
            if (i>0) //Not first
            {
                if (i<this->abscissae.size()-2) //Not last (abscissae got 1 new member, hence -2)
                {

                }
                else //Last
                {
                    auto im1=this->abscissae.size()-2;
                    auto wim1=this->weights.at(im1), wip1=this->weights.at(0), ai=0.0, bi=this->abscissae.at(0), ci=0.0;

                }
            }
            else //First
            {
                this->initialize_members();
            }
            */
        }

        template <class Generator>
        double operator()(Generator &g) const
        {
            bool bugged;
            double ret_value;
            do // while (!bugged)
            {
                bugged = false;
                try
                {
                    // std::unique_lock lock{this->m};
                    // rand is between 0 and normalization
                    double rand = static_cast<double>(g()) * this->normalization / (static_cast<double>(g.max() - g.min()));
                    // std::cout<<"here norm:"<<this->normalization<<" rand:"<<rand;
                    uint_fast16_t i = 0;
                    for (; i < this->cdf_weights.size() - 1; i++)
                    {
                        if (rand < this->cdf_weights.at(i + 1))
                        {
                            break;
                        }
                    }

                    // std::cout<<" i:"<<i;
                    auto ci = this->exp_factors.at(i);
                    auto cdfi = this->cdf_weights.at(i);
                    auto wi = this->weights.at(i);

                    // std::cout<<" ci:"<<ci<<" cdfi:"<<cdfi<<" wi:"<<wi<<std::endl;
                    ret_value = log(1 + ci * (rand - cdfi) / wi) / ci + this->abscissae.at(i);
                }
                catch (const std::exception &e)
                {
                    // std::cout << e.what() << " in piecewise_exp, trying again"<<std::endl;
                    bugged = true;
                }
            } while (bugged);

            return ret_value;
        }

    private:
        void initialize_members() noexcept
        {
            this->normalization = 0.0;
            auto wi = 0.0, wip1 = this->weights.at(0), ai = 0.0, bi = this->abscissae.at(0), ci = 0.0;
            this->exp_factors.clear();
            this->cdf_weights.clear();
            this->cdf_weights.emplace_back(0.0);
            for (uint_fast16_t ip1 = 1; ip1 < this->abscissae.size(); ip1++)
            {
                ai = bi;
                bi = this->abscissae.at(ip1);
                wi = wip1;
                wip1 = this->weights.at(ip1);

                ci = this->exp_factors.emplace_back(log(wip1 / wi) / (bi - ai));
                this->normalization += (wip1 - wi) / ci;
                this->cdf_weights.emplace_back(this->normalization);
            }
            // this->report();
        }
        void report() noexcept
        {
            std::cout << "normalization: " << this->normalization << std::endl
                      << "exp_factors: ";
            for (auto a : this->exp_factors)
            {
                std::cout << a << ' ';
            }
            std::cout << std::endl
                      << "cdf_weights: ";
            for (auto a : this->cdf_weights)
            {
                std::cout << a << ' ';
            }
            std::cout << std::endl
                      << "abscissae: ";
            for (auto a : this->abscissae)
            {
                std::cout << a << ' ';
            }
            std::cout << std::endl
                      << "weights: ";
            for (auto a : this->weights)
            {
                std::cout << a << ' ';
            }
            std::cout << std::endl
                      << std::endl;
        }

        double normalization{};
        std::vector<double> exp_factors{};
        std::vector<double> cdf_weights{};

        std::vector<double> abscissae{};
        std::vector<double> weights{};

        std::mutex m;
    };
    ars(const std::function<double(const double &)> &pdf_, const double &min, const double &max) noexcept;
    ars(const std::function<double(const double &)> &pdf_, const std::function<double(const double &)> &pdf_derivative_, const std::initializer_list<double> &init_abscissae, const std::initializer_list<double> &init_intersections, const std::initializer_list<double> &init_uk, const std::map<double, double> &init_log_values, const std::map<double, double> &init_log_der_values) noexcept;

    double throw_one_adaptive(std::mt19937 &random_generator);
    double throw_one_const(std::mt19937 &random_generator);

    double throw_one(std::mt19937 &random_generator)
    {
        if (this->adaptive)
        {
            return throw_one_adaptive(random_generator);
        }
        else
        {
            return throw_one_const(random_generator);
        }
    }
    std::vector<double> throw_n(const uint_fast32_t &n, std::mt19937 &random_generator) noexcept;
    bool is_adaptive() noexcept { return this->adaptive; }

private:
    void initialize_distribution() noexcept;
    void update_distribution() noexcept;

    void report() noexcept
    {
        std::cout << "log_values: ";
        for (auto a : this->log_values)
        {
            std::cout << '(' << a.first << ',' << a.second << ") ";
        }
        std::cout << std::endl
                  << "log_der_values: ";
        for (auto a : this->log_der_values)
        {
            std::cout << '(' << a.first << ',' << a.second << ") ";
        }
        std::cout << std::endl
                  << "abscissae: ";
        for (auto a : this->abscissae)
        {
            std::cout << a << ' ';
        }
        std::cout << std::endl
                  << "zk: ";
        for (auto a : this->zk)
        {
            std::cout << a << ' ';
        }
        std::cout << std::endl
                  << "uk: ";
        for (auto a : this->uk)
        {
            std::cout << a << ' ';
        }
        std::cout << std::endl
                  << std::endl;
    }

    bool adaptive{};

    const size_t absc_count_limit{}; // Maximum abscissae count

    const std::function<double(const double &)> log_pdf{};
    const std::function<double(const double &)> log_der{};

    std::map<double, double> log_values{};
    std::map<double, double> log_der_values{};
    std::vector<double> abscissae{};
    std::vector<double> zk{}; // Intersection points of the tangents of the logarithm of the pdf
    std::vector<double> uk{};

    std::uniform_real_distribution<double> unif_dist{};
    ars::piecewise_exponential_distribution sk{};

    std::mutex m;
};

#endif // ARS_HPP