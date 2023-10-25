//Copyright (c) 2023 Mikko Kuha

#ifndef TPP_HPP
#define TPP_HPP

#include <vector>

#include "nn_coll.hpp"
#include "nucleon.hpp"
#include "typedefs.hpp"

class Tpp_builder
{
public:
    Tpp_builder() : hotspots(false), proton_width_2(0.0), hs_width_2(0.0), func(nullptr) { }
    Tpp_builder(const double &proton_width_2_, const double &hotspot_width, const bool &hotspots_)
    : hotspots(hotspots_), proton_width_2(proton_width_2_)
    {
        if (hotspots)
        {
            this->hs_width_2 = std::pow(hotspot_width,2);
            this->func = std::function<double(const double&, const std::vector<double> *const)>
            ({[hs_width_2 = this->hs_width_2](const double &leave_zero = 0.0, const std::vector<double> *const hs_bs)
            {
                double Tpp = 0.0;
                uint_fast16_t n_hotspots2 = hs_bs->size(); //N*N
                for (auto bsquared : *hs_bs)
                {
                    Tpp += exp(-bsquared / (4 * hs_width_2)) / (40 * M_PI * hs_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
                }
                return Tpp/static_cast<double>(n_hotspots2);
            }});
        }
        else
        {
            this->hs_width_2 = 0.0;
            this->func = std::function<double(const double&, const std::vector<double> *const)>
            ({[proton_width_2 = this->proton_width_2](const double &bsquared, const std::vector<double> *const hs_bs = nullptr)
            {
                return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
            }});
        }
    }

    auto at
    (
        nn_coll &nn_pair
    ) noexcept -> double
    {
        if (hotspots)
        {
            auto vb2 = nn_pair.calculate_hs_bsquareds();
            return this->func(0.0, &vb2);
        }
        else
        {
            auto b2 = nn_pair.getcr_bsquared();
            return this->func(b2, nullptr);
        }
    }

    auto at
    (
        const double &b2
    ) noexcept -> double
    {
        if (hotspots)
        {
            //If at some point we need Tpp at singular b2, let's return the normal one proton Tpp 
            std::vector<double> vb2{b2 * this->hs_width_2 / this->proton_width_2};
            return this->func(0.0, &vb2) * this->hs_width_2 / this->proton_width_2;
        }
        else
        {
            return this->func(b2, nullptr);
        }
    }

    auto at
    (
        const std::vector<double> *const vb2
    ) noexcept -> double
    {
        if (!hotspots)
        {
            //If we got here, something is very wrong
            std::cout<<"HOTSPOT TPP CALLED WITH NO HOTSPOTS!!!!!!!"<<std::endl;
            return 0.0;
        }
        else
        {
            return this->func(0.0, vb2);
        }
    }

    auto calculate_sum_tpp
    (
        const nucleon &nuc, 
        const std::vector<nucleon> &nucleus
    ) noexcept -> double
    {
        double sum_tpp=0.0;
        if (!this->hotspots)
        {
            uint_fast16_t A = static_cast<uint_fast16_t>(nucleus.size());
            
            auto [x1, y1, z1] = nuc.co;
            
            std::vector<double> tpps(A, 0.0);
            
            for (uint_fast16_t i=0; i<A; i++)
            {
                if (nuc == nucleus.at(i)) // Do we calculate the effect of the same nucleon to itself?
                {
                    continue;
                }
                auto [x2, y2, z2] = nucleus.at(i).co;
                tpps.at(i) = this->func(pow(x1-x2,2) + pow(y1-y2,2), nullptr);
            }
            for (uint_fast16_t i=0; i<A; i++)
            {
                sum_tpp += tpps.at(i);
            }
        }
        else
        {
            uint_fast16_t A = static_cast<uint_fast16_t>(nucleus.size());
            uint_fast16_t N_hs = static_cast<uint_fast16_t>(nucleus.at(0).hotspots.size());
            
            auto hotspots = nuc.hotspots;
            
            std::vector<double> tpps(A, 0.0);
            
            for (uint_fast16_t i=0; i<A; i++)
            {
                if (nuc == nucleus.at(i)) // Do we calculate the effect of the same nucleon to itself?
                {
                    continue;
                }
                std::vector<double> b_squareds(N_hs*N_hs, 0.0);
                uint_fast16_t ind = 0;
                for (auto hs1 : hotspots)
                {
                    auto [x1, y1, z1] = hs1.co;
                    for (auto hs2 : nucleus.at(i).hotspots)
                    {
                        auto [x2, y2, z2] = hs2.co;
                        b_squareds.at(ind++) = pow(x1-x2,2) + pow(y1-y2,2);
                    }
                }

                tpps.at(i) = this->func(0.0, &b_squareds);
            }
            for (uint_fast16_t i=0; i<A; i++)
            {
                sum_tpp += tpps.at(i);
            }
        }
        return sum_tpp;
    }

    auto calculate_TA
    (
        const double &x,
        const double &y, 
        const std::vector<nucleon> &nucleus
    ) noexcept -> double
    {
        double sum_tp=0.0;
        if (!this->hotspots)
        {
            uint_fast16_t A = static_cast<uint_fast16_t>(nucleus.size());
            
            std::vector<double> tps(A, 0.0);
            
            for (uint_fast16_t i=0; i<A; i++)
            {
                auto [x2, y2, z2] = nucleus.at(i).co;
                // T^i_p(x) = 2*T_pp(2*(x-x_i)^2)
                tps.at(i) = 2.0 * this->func(2.0*(pow(x-x2,2) + pow(y-y2,2)), nullptr);
            }
            for (uint_fast16_t i=0; i<A; i++)
            {
                sum_tp += tps.at(i);
            }
        }
        else
        {
            uint_fast16_t A = static_cast<uint_fast16_t>(nucleus.size());
            uint_fast16_t N_hs = static_cast<uint_fast16_t>(nucleus.at(0).hotspots.size());
            
            std::vector<double> tps(A, 0.0);
            
            for (uint_fast16_t i=0; i<A; i++)
            {
                std::vector<double> b_squareds(N_hs, 0.0);
                uint_fast16_t ind = 0;
                for (auto hs : nucleus.at(i).hotspots)
                {
                    auto [x2, y2, z2] = hs.co;
                    b_squareds.at(ind++) = 2.0*(pow(x-x2,2) + pow(y-y2,2));
                }

                tps.at(i) = 2.0 * this->func(0.0, &b_squareds);
            }
            for (uint_fast16_t i=0; i<A; i++)
            {
                sum_tp += tps.at(i);
            }
        }
        return sum_tp;
    }

    auto calculate_T_AA_0
    (
        const nucleus_generator::nucleus_params &nuc_params,
        const double &rel_tolerance,
        std::shared_ptr<ars> radial_sampler,
        const bool verbose
    ) -> double
    {
        std::vector<uint_fast8_t> block_indexes(200);
        std::iota(block_indexes.begin(), block_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        
        uint_fast8_t block_amount = 50;
        std::vector<double> block_averages(block_amount);

        #pragma omp parallel
        {
            auto eng = std::make_shared<std::mt19937>(static_cast<ulong>((omp_get_thread_num() + 1))*static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));

            #pragma omp for
            for (uint_fast8_t block_index=0; block_index<block_amount; block_index++)
            {
                std::vector<double> T_AA_0s(200);
                for (auto it = block_indexes.begin(); it < block_indexes.end(); it++)
                {
                    std::vector<nucleon> nucl, other;
                    std::vector<uint_fast8_t> nucl_indexes(nuc_params.A);
                    std::iota(nucl_indexes.begin(), nucl_indexes.end(), 0);
                    std::vector<double> sum_tpps(nuc_params.A);

                    nucl = nucleus_generator::generate_nucleus
                        (
                            nuc_params,
                            false,
                            0.0, 
                            0.0, 
                            eng, 
                            radial_sampler
                        );

                    for (auto itt = nucl_indexes.begin(); itt < nucl_indexes.end(); itt++)
                    {
                        sum_tpps[*itt] = this->calculate_sum_tpp(nucl[*itt], nucl);
                    }

                    T_AA_0s[*it] = std::reduce(sum_tpps.begin(), sum_tpps.end(), 0.0);
                }

                block_averages[block_index] = std::reduce(T_AA_0s.begin(), T_AA_0s.end(), 0.0) / 200.0;
            }
        }

        double ave_T_AA_0 = std::reduce(block_averages.begin(), block_averages.end(), 0.0) / block_amount;

        std::cout<<"A-configs calculated: "<<block_amount*200<<std::endl;

        return ave_T_AA_0;
    }

    auto calculate_R_c_table
    (
        const nucleus_generator::nucleus_params &nuc_params,
        const double &rel_tolerance,
        std::shared_ptr<ars> radial_sampler,
        const bool verbose
    ) -> std::tuple<std::array<double, 25>, std::array<double, 25> >
    {
        std::vector<uint_fast8_t> block_indexes(200);
        std::iota(block_indexes.begin(), block_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        
        std::array<double, 25> c_vector{-15.0,-14.0,-13.0,-12.0,-11.0,-10.0,-9.0,-8.0,-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,10.0,15.0,20.0,25.0};
        std::array<double, 25> ave_R_vector{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

        uint_fast8_t block_amount = 50;
        std::array<std::array<double, 50>, 25> block_averages;

        #pragma omp parallel
        {
            auto eng = std::make_shared<std::mt19937>(static_cast<ulong>((omp_get_thread_num() + 1))*static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));

            #pragma omp for
            for (uint_fast8_t block_index=0; block_index<block_amount; block_index++)
            {
                std::array<std::array<double, 200>, 25> R_vectors;
                for (auto it = block_indexes.begin(); it < block_indexes.end(); it++)
                {
                    std::vector<nucleon> nucl, other;
                    std::vector<uint_fast8_t> nucl_indexes(nuc_params.A);
                    std::iota(nucl_indexes.begin(), nucl_indexes.end(), 0);
                    std::array<std::vector<double>, 25> sum_tpps;
                    for (uint_fast8_t i=0; i<25; i++)
                    {
                        sum_tpps[i] = std::vector<double>(nuc_params.A);
                    }

                    nucl = nucleus_generator::generate_nucleus
                        (
                            nuc_params,
                            false,
                            0.0, 
                            0.0, 
                            eng, 
                            radial_sampler
                        );

                    for (auto itt = nucl_indexes.begin(); itt < nucl_indexes.end(); itt++)
                    {
                        auto dummy = this->calculate_sum_tpp(nucl[*itt], nucl);
                        for (uint_fast8_t i=0; i<25; i++)
                        {
                            sum_tpps[i][*itt] = std::exp(c_vector[i]*dummy);
                        }
                    }

                    for (uint_fast8_t i=0; i<25; i++)
                    {
                        R_vectors[i][*it] = std::reduce(sum_tpps[i].begin(), sum_tpps[i].end(), 0.0) / static_cast<double>(nuc_params.A);
                    }
                }

                for (uint_fast8_t i=0; i<25; i++)
                {
                    block_averages[i][block_index] = std::reduce(R_vectors[i].begin(), R_vectors[i].end(), 0.0) / 200.0;
                }
            }
        }

        for (uint_fast8_t i=0; i<25; i++)
        {
            ave_R_vector[i] = std::reduce(block_averages[i].begin(), block_averages[i].end(), 0.0) / block_amount;
        }
        std::cout<<"A-configs calculated: "<<block_amount*200<<std::endl;

        return std::make_tuple(ave_R_vector,c_vector);
    }

private:
    bool hotspots{false};
    double proton_width_2{0.0};
    double hs_width_2{0.0};
    std::function<double(const double&, const std::vector<double> *const)> func{nullptr};
};

#endif // TPP_HPP