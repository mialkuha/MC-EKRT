//Copyright (c) 2023 Mikko Kuha

#ifndef PDF_HPP
#define PDF_HPP

#include <array>
#include <functional>
#include <omp.h>
#include <string>
#include <tuple>

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#include "eps09.h"
#include "LHAPDF/GridPDF.h"
#pragma GCC diagnostic pop

#include "linear_interpolator.hpp"
#include "nucleus_generator.hpp"
#include "Tpp_builder.hpp"
#include "typedefs.hpp"

class pdf_builder
{
public:
    pdf_builder
    (
        const std::string &p_pdf_name,
        const int &p_pdf_setnumber,
        const std::string &n_pdf_name,
        const int &n_pdf_setnumber,
        const double &snpdf_spatial_cutoff,
        const double &snpdf_tAA_0,
        const std::shared_ptr<Tpp_builder> Tpp,
        const nucleus_generator::nucleus_params &nuc_params,
        const bool &projectile_with_npdfs_,
        const bool &target_with_npdfs_,
        const bool &isoscalar_projectile_,
        const bool &isoscalar_target_,
        const bool &npdfs_spatial_,
        const bool &snPDFs_linear_,
        const bool &only_protons_,
        const int &npdf_setnumber_,
        const uint_fast16_t &A_,
        const uint_fast16_t &B_,
        const uint_fast16_t &ZA_,
        const uint_fast16_t &ZB_,
        const bool &verbose = false
    )
    : projectile_with_npdfs(projectile_with_npdfs_),
      target_with_npdfs(target_with_npdfs_),
      isoscalar_projectile(isoscalar_projectile_),
      isoscalar_target(isoscalar_target_),
      npdfs_spatial(npdfs_spatial_),
      snPDFs_linear(snPDFs_linear_),
      only_protons(only_protons_),
      npdf_setnumber(npdf_setnumber_),
      A(A_),
      B(B_),
      ZA(ZA_),
      ZB(ZB_),
      c_A_func(linear_interpolator(std::vector<double>{0.0},std::vector<double>{0.0}))
    {
        this->p_pdf = std::make_shared<LHAPDF::GridPDF>(p_pdf_name, p_pdf_setnumber);
        //this->n_pdf = std::make_shared<LHAPDF::GridPDF>(n_pdf_name, n_pdf_setnumber);

        if (this->npdfs_spatial && this->snPDFs_linear)
        {
            double tAA_0 = snpdf_tAA_0;
            if (tAA_0 == 0.0)
            {    
                std::cout<<"Calculating T_AA(0)"<<std::endl;

                tAA_0 = Tpp->calculate_T_AA_0
                (
                    nuc_params,
                    1e-5,
                    verbose
                );

                std::cout<<"Calculated T_AA(0) = "<< tAA_0 <<std::endl;
            }
            //r=(1+cT)
            //c=A*(R-1)/TAA(0)
            this->rA_spatial = [=, this](double const &r, double const &sum_tppa)
            {
                auto scaA = static_cast<double>(this->A) * sum_tppa / tAA_0;
                auto intA = 1.0 - scaA;
                auto dummy = r*scaA + intA;
                // return dummy; // no cutoff
                // return (dummy < 0.0 ) ? 0.0 : dummy; // cutoff at 1+cT < 0
                return (dummy < snpdf_spatial_cutoff ) ? snpdf_spatial_cutoff : dummy; // cutoff at 1+cT < spatial_cutoff
                // return (dummy < 1/static_cast<double>(NA) ) ? 1/static_cast<double>(NA) : dummy; // cutoff at 1+cT < A
                
            }; //r_s=1+c*sum(Tpp)

            this->rB_spatial = [=, this](double const &r, double const &sum_tppb)
            {
                auto scaB = static_cast<double>(this->B) * sum_tppb / tAA_0;
                auto intB = 1.0 - scaB;
                auto dummy = r*scaB + intB;
                // return dummy; // no cutoff
                // return (dummy < 0.0 ) ? 0.0 : dummy; // cutoff at 1+cT < 0
                return (dummy < snpdf_spatial_cutoff ) ? snpdf_spatial_cutoff : dummy; // cutoff at 1+cT < spatial_cutoff
                // return (dummy < 1/static_cast<double>(NB) ) ? 1/static_cast<double>(NB) : dummy; // cutoff at 1+cT < B
            };
        }
        else if (this->npdfs_spatial)
        {
            std::vector<double> Rs_as_vector, cs_as_vector;

            std::cout<<"Calculating R_A - c_A table"<<std::endl;
            auto [R_A_table, c_A_table] = Tpp->calculate_R_c_table
            (
                nuc_params,
                1e-5,
                verbose
            );
            std::cout<<"Done! {c,R} pairs:"<<std::endl;

            for (uint_fast8_t i=0; i<24; i++)
            {
                std::cout<<"{"<<c_A_table[i]<<","<<R_A_table[i]<<"},";
            }
            std::cout<<"{"<<c_A_table[24]<<","<<R_A_table[24]<<"}}"<<std::endl;
            
            for (uint_fast8_t i=0; i<25; i++)
            {
                Rs_as_vector.push_back(R_A_table[i]);
                cs_as_vector.push_back(c_A_table[i]);
            }

            this->c_A_func = linear_interpolator(Rs_as_vector, cs_as_vector);

            this->rA_spatial = [=, this](double const &r, double const &sum_tppa)
            {
                if (r>0.0)
                {
                    double c = this->c_A_func.value_at(r);
                    return std::exp(c * sum_tppa);
                }
                else
                {
                    return 0.0;
                }
            };

            this->rB_spatial = [=, this](double const &r, double const &sum_tppb)
            {
                if (r>0.0)
                {
                    double c = this->c_A_func.value_at(r);
                    return std::exp(c * sum_tppb);
                }
                else
                {
                    return 0.0;
                }
            };
        }
    }
    
    auto make_pdfs
    (
        const double &x1, 
        const double &x2, 
        const double &q2,
        const bool &target_neutron,
        const bool &projectile_neutron,
        const double &sum_tppa,
        const double &sum_tppb,
        const bool &average,
        const bool &max
    ) const noexcept -> std::tuple
        <
            std::array<double, 7>,
            std::array<double, 7>,
            std::array<double, 7>,
            std::array<double, 7>
        >
    {
        double rdv = 1.0, ruv = 1.0, rds = 1.0, rus = 1.0, rs = 1.0, rc = 1.0, rb = 1.0, rt = 1.0, rg = 1.0;
        double rd = 1.0, ru = 1.0;
        std::vector<double> xfx{13,0.0};
        this->p_pdf->xfxQ2(x1, q2, xfx);


        if (this->projectile_with_npdfs)
        {
            eps09(1, this->npdf_setnumber, static_cast<int>(this->A), x1, sqrt(q2), ruv, rdv, rus, rds, rs, rc, rb, rg);
            rd = rdv + (rds - rdv) * xfx[5] / xfx[7];
            ru = ruv + (rus - ruv) * xfx[4] / xfx[8];

            if (!average && this->npdfs_spatial)
            {
                rd = this->rA_spatial(rd, sum_tppa); ru = this->rA_spatial(ru, sum_tppa); rds = this->rA_spatial(rds, sum_tppa);
                rus = this->rA_spatial(rus, sum_tppa); rs = this->rA_spatial(rs, sum_tppa); rc = this->rA_spatial(rc, sum_tppa); 
                rb = this->rA_spatial(rb, sum_tppa); rt = this->rA_spatial(rt, sum_tppa); rg = this->rA_spatial(rg, sum_tppa);
            }

            if (max)
            {
                rd = (rd > 1.0)? rd: 1.0; ru = (ru > 1.0)? ru: 1.0; rds = (rds > 1.0)? rds: 1.0;
                rus = (rus > 1.0)? rus: 1.0; rs = (rs > 1.0)? rs: 1.0; rc = (rc > 1.0)? rc: 1.0; 
                rb = (rb > 1.0)? rb: 1.0; rt = (rt > 1.0)? rt: 1.0; rg = (rg > 1.0)? rg: 1.0;
            }
        }
        std::array<double, 7> f_i_x1 = 
        {
            rg * xfx[6],
            rd * xfx[7],
            ru * xfx[8],
            rs * xfx[9],
            rc * xfx[10],
            rb * xfx[11],
            rt * xfx[12]
        };
        std::array<double, 7> f_ai_x1 = 
        {
            rg * xfx[6],
            rds * xfx[5],
            rus * xfx[4],
            rs * xfx[3],
            rc * xfx[2],
            rb * xfx[1],
            rt * xfx[0]
        };


        this->p_pdf->xfxQ2(x2, q2, xfx);

        if (this->target_with_npdfs)
        {
            eps09(1, this->npdf_setnumber, static_cast<int>(this->B), x2, sqrt(q2), ruv, rdv, rus, rds, rs, rc, rb, rg);
            rd = rdv + (rds - rdv) * xfx[5] / xfx[7];
            ru = ruv + (rus - ruv) * xfx[4] / xfx[8];

            if (!average && this->npdfs_spatial)
            {
                rd = this->rB_spatial(rd, sum_tppb); ru = this->rB_spatial(ru, sum_tppb); rds = this->rB_spatial(rds, sum_tppb);
                rus = this->rB_spatial(rus, sum_tppb); rs = this->rB_spatial(rs, sum_tppb); rc = this->rB_spatial(rc, sum_tppb); 
                rb = this->rB_spatial(rb, sum_tppb); rt = this->rB_spatial(rt, sum_tppb); rg = this->rB_spatial(rg, sum_tppb);
            }

            if (max)
            {
                rd = (rd > 1.0)? rd: 1.0; ru = (ru > 1.0)? ru: 1.0; rds = (rds > 1.0)? rds: 1.0;
                rus = (rus > 1.0)? rus: 1.0; rs = (rs > 1.0)? rs: 1.0; rc = (rc > 1.0)? rc: 1.0; 
                rb = (rb > 1.0)? rb: 1.0; rt = (rt > 1.0)? rt: 1.0; rg = (rg > 1.0)? rg: 1.0;
            }
        }
        else
        {
            rd = 1.0; ru = 1.0; rdv = 1.0; ruv = 1.0; rds = 1.0; rus = 1.0; 
            rs = 1.0; rc = 1.0; rb = 1.0; rt = 1.0; rg = 1.0; 
        }
        std::array<double, 7> f_i_x2 = 
        {
            rg * xfx[6],
            rd * xfx[7],
            ru * xfx[8],
            rs * xfx[9],
            rc * xfx[10],
            rb * xfx[11],
            rt * xfx[12]
        };
        std::array<double, 7> f_ai_x2 = 
        {
            rg * xfx[6],
            rds * xfx[5],
            rus * xfx[4],
            rs * xfx[3],
            rc * xfx[2],
            rb * xfx[1],
            rt * xfx[0]
        };

        if (!this->only_protons && projectile_neutron)
        {
            double u, ub, d, db;
            d  = f_i_x1[1];
            u  = f_i_x1[2];
            db = f_ai_x1[1];
            ub = f_ai_x1[2];
            f_i_x1[1] = u;
            f_i_x1[2] = d;
            f_ai_x1[1] = ub;
            f_ai_x1[2] = db;
        }    
        if (!this->only_protons && target_neutron)
        {
            double u, ub, d, db;
            d  = f_i_x2[1];
            u  = f_i_x2[2];
            db = f_ai_x2[1];
            ub = f_ai_x2[2];
            f_i_x2[1] = u;
            f_i_x2[2] = d;
            f_ai_x2[1] = ub;
            f_ai_x2[2] = db;
        }

        return std::make_tuple(f_i_x1, f_i_x2, f_ai_x1, f_ai_x2);
    }

    auto alphasQ2(const double &q2) const -> const double { return this->p_pdf->alphasQ2(q2); };
    auto alphasQ(const double &q) const -> const double { return this->p_pdf->alphasQ(q); };
    auto num_flavors() const -> const uint_fast8_t { return static_cast<uint_fast8_t>(std::stoi(this->p_pdf->info().get_entry("NumFlavors"))); };
    auto set_index() const -> const std::string { return this->p_pdf->info().get_entry("SetIndex"); };
private:
    bool projectile_with_npdfs{true};
    bool target_with_npdfs{true};
    bool isoscalar_projectile{false};
    bool isoscalar_target{false};
    bool npdfs_spatial{true};
    bool snPDFs_linear{false};   
    bool only_protons{false};
    int npdf_setnumber{1};
    uint_fast16_t A{208};
    uint_fast16_t B{208};
    uint_fast16_t ZA{82};
    uint_fast16_t ZB{82};
    std::shared_ptr<LHAPDF::GridPDF> p_pdf{};
    std::shared_ptr<LHAPDF::GridPDF> n_pdf{};
    std::function<double (const double &, const double &)> rA_spatial{};
    std::function<double (const double &, const double &)> rB_spatial{};
    linear_interpolator c_A_func;

};

#endif // PDF_HPP