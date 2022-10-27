//Copyright (c) 2022 Mikko Kuha 
//TODO implement M-tree (https://github.com/erdavila/M-Tree ? Boost::Graph?)

#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <execution>
#include <fstream>
#include <gsl/gsl_sf_expint.h>
#include <iostream>
#include <mutex>
#include <random>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "ars.hpp"
#include "generic_helpers.hpp"
#include "high_level_calcs.hpp"
#include "histo.hpp"
#include "io_helpers.hpp"
#include "LHAPDF/GridPDF.h"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

auto find_sigma_jet_cutoff_p0
(
    momentum &kt02, 
    const momentum &mand_s, 
    const xsectval &target, 
    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    const pqcd::sigma_jet_params &jet_params, 
    const bool &verbose=true
) noexcept -> momentum
{
    xsectval sigma_jet=0.0;
    kt02 = 2.0;

    auto difference_to_target = [&](const momentum &_kt02)
    {
        return pqcd::calculate_sigma_jet(p_p_pdf, &mand_s, &_kt02, jet_params) - target;
    };

    helpers::secant_method(&kt02, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<sqrt(kt02)<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return kt02;
}

auto find_sigma_jet_cutoff_Q
(
    const momentum &kt02, 
    const momentum &mand_s, 
    const xsectval &target, 
    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    pqcd::sigma_jet_params jet_params, 
    const bool &verbose=true
) noexcept -> void
{
    xsectval sigma_jet=0.0;
    momentum scalar = 1.0;

    auto difference_to_target = [&](const momentum &scalar_)
    {
        momentum kt02dummy = 4.0;
        auto jet_params_ = pqcd::sigma_jet_params(
        /*d_params=                 */jet_params.d_params,
        /*scale_choice=             */jet_params.scale_c,
        /*scalar=                   */scalar_,
        /*use_ses=                  */jet_params.use_ses);
        auto kt02_ = find_sigma_jet_cutoff_p0(kt02dummy, mand_s, target, p_p_pdf, jet_params_);
        return kt02_ - kt02;
    };

    helpers::secant_method(&scalar, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<scalar<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return;
}

auto check_and_place_circle_among_others_with_overlaps
(
    std::tuple<double, double, double, uint16_t> cand_circle, 
    std::vector<std::tuple<double, double, double, uint16_t> > &final_circles, 
    const double &maximum_overlap,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> bool
{
    auto & [cand_x, cand_y, cand_r, cand_overlap] = cand_circle;
    std::vector<uint16_t*> overlaps_with;
    
    for (auto & [circ_x, circ_y, circ_r, circ_overlap] : final_circles)
    {
        if ( pow((circ_x-cand_x),2) + pow((circ_y-cand_y),2) < pow((circ_r+cand_r),2) )
        {
            cand_overlap += 1;
            if (cand_overlap > maximum_overlap)
            {
                if
                (
                    ((cand_overlap - maximum_overlap) >= 1) || 
                    ((cand_overlap - maximum_overlap) > unirand(*random_generator))
                )
                {
                    return false;
                }
            }
            if ((circ_overlap + 1) > maximum_overlap)
            {
                if 
                (
                    (circ_overlap >= maximum_overlap) || 
                    ((circ_overlap + 1 - maximum_overlap) > unirand(*random_generator))
                )
                {
                    return false;
                }
            }
            overlaps_with.push_back(&circ_overlap);
        }
    }
    
    for (auto & c : overlaps_with)
    {
        (*c)++;
    }

    final_circles.emplace_back(std::move(cand_circle));
        
    return true;
}

auto check_and_place_circle_among_others_with_disks
(
    std::tuple<double, double, double, uint16_t> cand_circle, 
    std::vector<std::tuple<double, double, double, uint16_t> > &final_circles
) noexcept -> bool
{
    auto & [cand_x, cand_y, cand_r, cand_overlap] = cand_circle;
    std::vector<uint16_t*> overlaps_with;
    
    for (auto & [circ_x, circ_y, circ_r, circ_overlap] : final_circles)
    {
        if ( pow((circ_x-cand_x),2) + pow((circ_y-cand_y),2) < pow((circ_r+cand_r),2) )
        {
            return false;
        }
    }

    final_circles.emplace_back(std::move(cand_circle));
        
    return true;
}

struct dijet_with_ns
{
    dijet_specs dijet;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    double tau;
    dijet_with_ns(dijet_specs dijet_, nucleon * pro_nucleon_, nucleon * tar_nucleon_)
        : dijet(std::move(dijet_)), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_)
    { tau = std::cosh(std::max(std::abs(dijet.y1),std::abs(dijet.y2)))/(dijet.kt*(1+std::exp(-std::abs(dijet.y1-dijet.y2))));}
};

auto throw_location_for_dijet
(
    const dijet_with_ns &cand,
    const spatial &proton_width,
    std::normal_distribution<double> &normal_dist, 
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> std::tuple<double, double, double>
{
    auto param = std::normal_distribution<double>::param_type{0., proton_width};
    auto dx = normal_dist(*random_generator,param);
    auto dy = normal_dist(*random_generator,param);
    auto dz = 0.0;

    auto retx = 0.5*(cand.pro_nucleon->co.x + cand.tar_nucleon->co.x + M_SQRT2*dx);
    auto rety = 0.5*(cand.pro_nucleon->co.y + cand.tar_nucleon->co.y + M_SQRT2*dy);

    return std::make_tuple(retx, rety, dz);
}

//Checks whether the saturation criterium allows the candidate to be added TODO
auto filter_collisions_saturation
(
    std::vector<dijet_with_ns> &candidates, 
    std::vector<dijet_with_coords> &final_candidates,
    const double &maximum_overlap,
    const spatial &proton_width,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator,
    const bool with_overlaps,
    const bool include_tata = false,
    const std::vector<nucleon> &pro = {}, 
    const std::vector<nucleon> &tar = {}, 
    const std::function<spatial(const spatial&)> Tpp = nullptr
) noexcept -> void
{
    std::vector<std::tuple<double, double, double, uint16_t> > final_circles;
    final_circles.reserve(candidates.size());
    final_candidates.reserve(final_candidates.size()+candidates.size());

    std::normal_distribution<double> normal_dist(0,0);

    //std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the biggest kt is first
    //    [](dijet_with_ns &s1, dijet_with_ns &s2) { return (s1.dijet.kt > s2.dijet.kt); });
    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the smallest tau is first
        [](dijet_with_ns &s1, dijet_with_ns &s2) { return (s1.tau < s2.tau); });

    if (!with_overlaps)
    {
        for (auto & cand : candidates)
        {
            auto [cand_x, cand_y, cand_z] = throw_location_for_dijet(cand, proton_width, normal_dist, random_generator);
            auto cand_circle = std::make_tuple(cand_x, cand_y, 1/(cand.dijet.kt*maximum_overlap*FMGEV), 1);
            if (check_and_place_circle_among_others_with_disks(std::move(cand_circle), final_circles))
            {
                double tata = 0.0;
                if (include_tata)
                {
                    nucleon dummy{coords{cand_x, cand_y, cand_z}, 0};
                    tata = calcs::calculate_sum_tpp(dummy, pro, Tpp) * calcs::calculate_sum_tpp(dummy, tar, Tpp);
                }
                final_candidates.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.tau, tata});
            }
        }
    }
    else
    {
        for (auto & cand : candidates)
        {
            auto [cand_x, cand_y, cand_z] = throw_location_for_dijet(cand, proton_width, normal_dist, random_generator);
            auto cand_circle = std::make_tuple(cand_x, cand_y, 1/cand.dijet.kt, 1);
            if (check_and_place_circle_among_others_with_overlaps(std::move(cand_circle), final_circles, maximum_overlap, unirand, random_generator))
            {
                double tata = 0.0;
                if (include_tata)
                {
                    nucleon dummy{coords{cand_x, cand_y, cand_z}, 0};
                    tata = calcs::calculate_sum_tpp(dummy, pro, Tpp) * calcs::calculate_sum_tpp(dummy, tar, Tpp);
                }
                final_candidates.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.tau, tata});
            }
        }
    }
    
    final_candidates.shrink_to_fit();
}

struct colls_with_ns
{
    momentum kt;
    rapidity y1;
    rapidity y2;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    dijet_specs * dijet;
    double tau;
    colls_with_ns(momentum kt_, rapidity y1_, rapidity y2_, nucleon * pro_nucleon_, nucleon * tar_nucleon_, dijet_specs * dijet_)
        : kt(kt_), y1(y1_), y2(y2_), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_), dijet(dijet_)
    { tau = std::cosh(std::max(std::abs(y1),std::abs(y2)))/(kt*(1+std::exp(-std::abs(y1-y2))));}
};
//Empties all of the binary_collisions
auto filter_collisions_MC 
(
    std::vector<nn_coll> &binary_collisions,
    std::vector<dijet_with_ns> &final_candidates,
    const momentum sqrt_s,
    const bool deplete_nucleons
) noexcept -> void
{
    std::unordered_map<nucleon*, double> x1s;
    std::unordered_map<nucleon*, double> x2s;
    std::unordered_set<nucleon*> depleted_pro;
    std::unordered_set<nucleon*> depleted_tar;

    std::vector<colls_with_ns> collision_candidates;
    collision_candidates.reserve(binary_collisions.size()*10);

    for (auto &col : binary_collisions)
    {
        for (auto &dij : col.dijets)
        {
            collision_candidates.emplace_back(dij.kt, dij.y1, dij.y2, col.projectile, col.target, &dij);
        }
    }
    collision_candidates.shrink_to_fit();

    //std::sort(collision_candidates.begin(), collision_candidates.end(), //Sort the candidates so that the one with the biggest kt is first
    //          [](colls_with_ns &s1, colls_with_ns &s2) { return (s1.kt > s2.kt); });
    std::sort(collision_candidates.begin(), collision_candidates.end(), //Sort the candidate events so that the one with the smallest tau is first
              [](colls_with_ns &s1, colls_with_ns &s2) { return (s1.tau < s2.tau); });

    for (auto & cand : collision_candidates)
    {
        if (deplete_nucleons && (depleted_pro.contains(cand.pro_nucleon) || depleted_tar.contains(cand.tar_nucleon)))
        {
            continue;
        }

        bool discard = false;
        auto x1 = (cand.kt / sqrt_s) * (exp(cand.y1) + exp(cand.y2));
        auto x2 = (cand.kt / sqrt_s) * (exp(-cand.y1) + exp(-cand.y2));

        auto i_x1_sum = x1s.find(cand.pro_nucleon);
        if (i_x1_sum != x1s.end())
        {
            i_x1_sum->second += x1;

            if (i_x1_sum->second > 1.0)
            {
                discard = true;
                if (deplete_nucleons)
                {
                    depleted_pro.insert(cand.pro_nucleon);
                }
            }
        }
        else
        {
            x1s.insert({cand.pro_nucleon, x1});
        }

        auto i_x2_sum = x2s.find(cand.tar_nucleon);
        if (i_x2_sum != x2s.end())
        {
            i_x2_sum->second += x2;

            if (i_x2_sum->second > 1.0)
            {
                discard = true;
                if (deplete_nucleons)
                {
                    depleted_tar.insert(cand.tar_nucleon);
                }
            }
        }
        else
        {
            x2s.insert({cand.tar_nucleon, x2});
        }

        if (!discard)
        {
            final_candidates.emplace_back(std::move(*(cand.dijet)), cand.pro_nucleon, cand.tar_nucleon);
        }
    }

    //std::cout << "candidates filtered: " << collision_candidates.size() - final_candidates.size() << " out of " << collision_candidates.size() << std::endl;

    binary_collisions.clear();
    //for (auto &col : binary_collisions)
    //{
    //    col.dijets.erase(col.dijets.begin(), col.dijets.end());
    //}
}


//Empties all of the binary_collisions TODO
auto filter_collisions_local_MC 
(
    std::vector<nn_coll> &binary_collisions,
    std::vector<dijet_with_ns> &final_candidates,
    const momentum &sqrt_s,
    const xsectval &sigma_inel
) noexcept -> void
{
    /*std::unordered_map<nucleon*, double> x1s;
    std::unordered_map<nucleon*, double> x2s;
    std::unordered_set<nucleon*> depleted_pro;
    std::unordered_set<nucleon*> depleted_tar;

    std::vector<colls_with_ns> collision_candidates;
    collision_candidates.reserve(binary_collisions.size()*10);

    for (auto &col : binary_collisions)
    {
        for (auto &dij : col.dijets)
        {
            collision_candidates.emplace_back(dij.kt, dij.y1, dij.y2, col.projectile, col.target, &dij);
        }
    }
    collision_candidates.shrink_to_fit();

    std::sort(collision_candidates.begin(), collision_candidates.end(), //Sort the candidates so that the one with the biggest kt is first
              [](colls_with_ns &s1, colls_with_ns &s2) { return (s1.kt > s2.kt); });

    for (auto & cand : collision_candidates)
    {
        if (deplete_nucleons && (depleted_pro.contains(cand.pro_nucleon) || depleted_tar.contains(cand.tar_nucleon)))
        {
            continue;
        }

        bool discard = false;
        auto x1 = (cand.kt / sqrt_s) * (exp(cand.y1) + exp(cand.y2));
        auto x2 = (cand.kt / sqrt_s) * (exp(-cand.y1) + exp(-cand.y2));

        auto i_x1_sum = x1s.find(cand.pro_nucleon);
        if (i_x1_sum != x1s.end())
        {
            i_x1_sum->second += x1;

            if (i_x1_sum->second > 1.0)
            {
                discard = true;
                if (deplete_nucleons)
                {
                    depleted_pro.insert(cand.pro_nucleon);
                }
            }
        }
        else
        {
            x1s.insert({cand.pro_nucleon, x1});
        }

        auto i_x2_sum = x2s.find(cand.tar_nucleon);
        if (i_x2_sum != x2s.end())
        {
            i_x2_sum->second += x2;

            if (i_x2_sum->second > 1.0)
            {
                discard = true;
                if (deplete_nucleons)
                {
                    depleted_tar.insert(cand.tar_nucleon);
                }
            }
        }
        else
        {
            x2s.insert({cand.tar_nucleon, x2});
        }

        if (!discard)
        {
            final_candidates.emplace_back(std::move(*(cand.dijet)), cand.pro_nucleon, cand.tar_nucleon);
        }
    }

    //std::cout << "candidates filtered: " << collision_candidates.size() - final_candidates.size() << " out of " << collision_candidates.size() << std::endl;

    binary_collisions.clear();*/
    //for (auto &col : binary_collisions)
    //{
    //    col.dijets.erase(col.dijets.begin(), col.dijets.end());
    //}
}

auto filter_end_state
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_with_coords> &filtered_scatterings,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator,
    const bool mom_cons = false,
    const bool local_mom_cons = false,
    const bool saturation = false,
    const bool deplete_nucleons = false,
    const bool saturation_with_overlap = false,
    const momentum sqrt_s = 0,
    const double &maximum_overlap = 2.0,
    const spatial &proton_width = 1.0,
    const xsectval &sigma_inel = 1.0,
    const bool include_tata = false,
    const std::vector<nucleon> &pro = {}, 
    const std::vector<nucleon> &tar = {}, 
    const std::function<spatial(const spatial&)> Tpp = nullptr
) noexcept -> void
{
    std::vector<dijet_with_ns> candidates;
    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess

    std::vector<dijet_specs> placeholder_filtered;
    
    if (mom_cons) 
    {
        if (local_mom_cons)
        {
            filter_collisions_local_MC(binary_collisions, candidates, sqrt_s, sigma_inel);
        }
        else
        {
            filter_collisions_MC(binary_collisions, candidates, sqrt_s, deplete_nucleons);
        }
    }
    else
    {
        for (auto &col : binary_collisions) //all of the end states of the binary events are now candidate events
        {
            for (auto & dij : col.dijets)
            {
                candidates.emplace_back(std::move(dij), col.projectile, col.target);
            }
        }
        binary_collisions.clear();
    }
    candidates.shrink_to_fit();
    

    if (saturation)
    {
        filter_collisions_saturation(candidates, filtered_scatterings, maximum_overlap, proton_width, unirand, random_generator, saturation_with_overlap, include_tata, pro, tar, Tpp);
    }
    else
    {
        for (auto &cand : candidates) //all of the candidates are passed to filtered
        {
            std::normal_distribution<double> normal_dist(0,0);
            double tata = 0.0;
            auto [cand_x, cand_y, cand_z] = throw_location_for_dijet(cand, proton_width, normal_dist, random_generator);
            if (include_tata)
            {
                nucleon dummy{coords{cand_x, cand_y, cand_z}, 0};
                tata = calcs::calculate_sum_tpp(dummy, pro, Tpp) * calcs::calculate_sum_tpp(dummy, tar, Tpp);
            }
            filtered_scatterings.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.tau, tata});
            //filtered_scatterings.push_back(candidate.dijet);
        }
    }
    
    //std::cout << "filtered_scatterings.size() =  " << filtered_scatterings.size() << std::endl;
}

volatile std::atomic_bool user_aborted = false;
volatile std::atomic_bool g_bug_bool = false;
void abort_handler(int num)
{
    user_aborted = true;
}
void term_handler(int num)
{
    std::cout<<std::endl<<"Terminate called "<<num<<std::endl;
}

int main(int argc, char** argv)
{ 
    /*bool g_is_pp;
    bool g_is_pa;
    bool g_is_aa;
    bool g_use_npdfs;
    bool g_use_snpdfs;
    bool g_is_mom_cons;
    bool g_is_mom_cons_new;
    bool g_are_ns_depleted;
    bool g_is_saturation;
    std::string name_postfix;*/

    auto 
    [
        name_postfix,
        g_is_pp,
        g_is_pa,
        g_is_aa,
        g_use_npdfs,
        g_use_snpdfs,
        g_is_mom_cons,
        g_is_mom_cons_new,
        g_are_ns_depleted,
        g_is_saturation,
        g_is_sat_overlap,
        g_is_mc_glauber,
        desired_N_events,
        b_max,
        pt0,
        K_factor,
        K_sat
    ] = io::read_conf(std::string(argv[1]));

    /*switch (argv[1][0])
    {
    case '0':
        name_postfix = "_pA_2500k_mb_PDF";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '1':
        name_postfix = "_pA_2500k_mb_nPDF";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '2':
        name_postfix = "_pA_2500k_mb_snPDF";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = true;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '3':
        name_postfix = "_pA_2500k_mb_PDF_MC";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '4':
        name_postfix = "_pA_2500k_mb_nPDF_MC";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '5':
        name_postfix = "_pA_2500k_mb_snPDF_MC";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = true;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '6':
        name_postfix = "_pA_2500k_mb_PDF_MC_ND";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = true;
        g_is_saturation = false;
        break;
    case '7':
        name_postfix = "_pA_2500k_mb_nPDF_MC_ND";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = true;
        g_is_saturation = false;
        break;
    case '8':
        name_postfix = "_pA_2500k_mb_snPDF_MC_ND";
        g_is_pp = false;
        g_is_pa = true;
        g_is_aa = false;
        g_use_npdfs = true;
        g_use_snpdfs = true;
        g_is_mom_cons = false;
        g_is_mom_cons_new = true;
        g_are_ns_depleted = true;
        g_is_saturation = false;
        break;
    default:
        break;
    }*/

    //A lot of printing
    bool verbose = false;

    if (verbose) std::cout<<"Initializing..."<<std::flush;
    
    //General parameters for the simulation
    const bool    read_nuclei_from_file    = false, 
                  read_events_from_file    = false,
                  read_sigmajets_from_file = true,
                  end_state_filtering      = true, 
                  varying_inel             = true, 
                  save_events              = false, 
                  include_tata             = true/*, 
                  average_spatial_taas     = false*/;
    std::string   event_file_name = "event_log_"+name_postfix+".dat";
    uint32_t      AA_events             = 0;
    const spatial b_min                 = 0;
    std::mutex AA_events_mutex; 
    std::cout<<std::endl<<std::endl<<"Doing the run "<<name_postfix<<std::endl;
    //auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(1));
    auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand{0.0, 1.0};

    //Parameters for the nuclei
    //uint8_t samplersN = 32;
    const spatial rad_min=0,
                  rad_max=30;
    const std::function<double(const double&)> rad_pdf{[](const double & x)
    {
        return x*x/(1+exp((x-6.624)/0.549));
    }};
    auto radial_sampler{std::make_shared<ars>(rad_pdf, rad_min, rad_max)};
    do //while (radial_sampler->is_adaptive())
    {
        radial_sampler->throw_one(*eng);
    } while (radial_sampler->is_adaptive());
    //std::vector<std::shared_ptr<ars> > sampler_pointers;
    //std::mutex radial_sampler_mutex;
    //for (uint8_t i=0; i<samplersN; i++)
    //{
    //    sampler_pointers.push_back(std::make_shared<ars>(rad_pdf, rad_min, rad_max));
    //}

    nucleus_generator::nucleus_params nuc_params = 
    {
        /* .NA=                   */(g_is_aa)? 208u : 1u, //Pb 
        /* .ZA=                   */(g_is_aa)? 82u : 1u, //Pb 
        /* .NB=                   */(g_is_pp)? 1u : 208u, //Pb 
        /* .ZB=                   */(g_is_pp)? 1u : 82u, //Pb 
        /* .min_distance=         */0.4, 
        /* .shift_cms=            */true, 
        /* .correct_overlap_bias= */true
    };
    
    //Parameters for the hard collisions
    const spatial proton_width = 0.573;
    const spatial proton_width_2 = pow(0.573, 2);
    const std::function<spatial(const spatial&)> Tpp{[&proton_width_2](const spatial &bsquared)
    {
        return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
    }}; 
    xsectval sigma_inel_for_glauber       = 70;//mb
    const momentum sqrt_s                 = 5020;//GeV
    const momentum mand_s                 = pow(sqrt_s, 2);//GeV^2
    momentum kt0                          = pt0;//2.728321;//GeV
    momentum kt02                         = pow(kt0, 2);//GeV^2
    rapidity ylim                         = static_cast<rapidity>(log(sqrt_s / kt0));
    auto p_pdf = std::make_shared<LHAPDF::GridPDF>("CT14lo", 0);

    std::array<std::vector<io::Coll>, 6> collisions_for_reporting;
    std::array<std::vector<io::Coll>, 6> collisions_for_reporting_midrap;
    std::mutex colls_mutex; 
    std::vector<double> kt_bins{helpers::loglinspace(kt0, 100, 16)};
    auto kt_v_dummy = helpers::loglinspace(200, sqrt_s/2.0, 5);
    for (auto ktdummy : kt_v_dummy)
    {
        kt_bins.push_back(ktdummy);
    }
    const std::vector<double> y_bins{helpers::linspace(-ylim, ylim, 40)};
    const std::vector<double> b_bins{helpers::linspace(b_min, b_max, 21)};
    std::vector<double> et_bins{helpers::loglinspace(2*kt0, 30000, 31)};

    //sigma_jet parameters
    /*const bool read_sigmajets_from_file = false;*/
    pqcd::diff_sigma::params diff_params = pqcd::diff_sigma::params(
    /*projectile_with_npdfs=    */(g_is_aa && g_use_npdfs),
    /*target_with_npdfs=        */(!g_is_pp && g_use_npdfs),
    /*isoscalar_projectile=     */false,
    /*isoscalar_target=         */false,
    /*npdfs_spatial=            */g_use_snpdfs,
    /*npdf_setnumber=           */1,
    /*K_factor=                 */K_factor,
    /*A=                        */(g_is_aa)? 208 : 1, //Pb 
    /*B=                        */(g_is_pp)? 1 : 208, //Pb
    /*ZA=                       */(g_is_aa)? 82 : 1,  //Pb
    /*ZB=                       */(g_is_pp)? 1 : 82   //Pb
    /*p_n_pdf=                  */
    /*rA_spatial=               */
    /*rB_spatial=               */);
    pqcd::sigma_jet_params jet_params = pqcd::sigma_jet_params(
    /*d_params=                 */diff_params,
    /*scale_choice=             */pqcd::scaled_from_kt,
    /*scalar=                   */1.0,
    /*use_ses=                  */false);

    //calculate_sigma_1jet_analytical(mand_s, kt_bins, y_bins, p_pdf, &jet_params, name_postfix);
    //calculate_sigma_jet_analytical(mand_s, kt_bins, y_bins, p_pdf, &jet_params);
    //calculate_sigma_dijet_analytical(mand_s, kt_bins, y_bins, p_pdf, &jet_params, name_postfix);
//std::vector<double> mands = {2., 5., 10., 20., 50.,100.,200.,500.,1000.,2000.,5000.,};
//for (auto & m : mands) m = m*m;
//std::vector<double> kt02ss = {1., 1., 1., 1., 1.,1.,1.,1.,1.,1.,1.,};
//std::vector<double> inels = {1., 1., 1., 1., 1.,1.,1.,1.,1.,1.,1.,};
//for (int i=0; i<11 ; i++) find_sigma_jet_cutoff(kt02ss[i], mands[i], 124.6635, p_pdf, jet_params, true);
//std::cout<<kt02<<std::endl;
//xsectval s_jet = 124.6635;
//xsectval target_s_jet;
//auto scalars = {1.0, 1/2.0, 1/3.0, 1/4.0};
//auto kt02s = {pow(2.11167, 2), pow(3.45886, 2), pow(4.33213, 2), pow(5.38004, 2), pow(6.63009, 2)};
//std::vector<double> kt02s{helpers::loglinspace(0.25, 100.0, 150)};
//std::ofstream out_file;
//out_file.open("p02_sjet.csv", std::ios::out);
//
//find_sigma_jet_cutoff_p0(kt02, mand_s, target_s_jet, p_pdf, jet_params);
//
//for (auto s : scalars)
//{
//    target_s_jet = s*s_jet;
//    find_sigma_jet_cutoff_p0(kt02, mand_s, target_s_jet, p_pdf, jet_params);
//    kt02s.push_back(find_sigma_jet_cutoff_p0(kt02, mand_s, target_s_jet, p_pdf, jet_params));
//    std::cout<<std::endl;
//}
//
//for (auto k : kt02s)
//{
//    auto res = pqcd::calculate_sigma_jet(p_pdf, &mand_s, &k, jet_params);
//    std::cout<<k<<','<<res<<std::endl;
//    out_file<<k<<','<<res<<std::endl;
//}
//for (auto k : kt02s)
//{
//    find_sigma_jet_cutoff_Q(k, mand_s, s_jet, p_pdf, jet_params);
//    std::cout<<std::endl;
//}
//
//return 0;

    auto
    [
        dijet_norm,
        power_law,
        sigma_jet,
        mand_s_array,
        envelope_maximum,
        sqrt_s_array
    ] = 
    calcs::prepare_sigma_jets
    (
        (diff_params.projectile_with_npdfs || diff_params.target_with_npdfs),
        g_use_snpdfs,
        g_is_mom_cons,
        read_sigmajets_from_file,
        p_pdf, 
        mand_s,
        sqrt_s,
        kt02, 
        kt0,
        jet_params
    );
    
    if (varying_inel)
    {
        auto dummy = std::get<xsectval>(sigma_jet) / (4 * M_PI * proton_width_2);
        dummy = dummy / 10; //mb -> fm²

        sigma_inel_for_glauber = (4 * M_PI * proton_width_2) * (M_EULER + std::log(dummy) + gsl_sf_expint_E1(dummy)) * 10; //1 mb = 0.1 fm^2
        std::cout<<std::endl<<sigma_inel_for_glauber<<std::endl;
    }

    AA_collision_params coll_params
    {
    /*mc_glauber_mode=          */g_is_mc_glauber,
    /*pp_scattering=            */g_is_pp,
    /*pA_scattering=            */g_is_pa,
    /*spatial_pdfs=             */g_use_snpdfs,
    /*spatial_averaging=        */false,
    /*calculate_end_state=      */true,
    /*reduce_nucleon_energies=  */g_is_mom_cons,
    /*sigma_inel_for_glauber=   */sigma_inel_for_glauber,
    /*Tpp=                      */Tpp,
    /*normalize_to=             */B2_normalization_mode::inelastic,
    /*sqrt_s=                   */sqrt_s,
    /*energy_threshold=         */kt0
    };

//std::cout<<std::get<xsectval>(sigma_jet)<<std::endl;
//return 0;

    //find_sigma_jet_cutoff(kt02, mand_s, 124.6635, p_pdf, jet_params, true);
    //std::cout<<kt02<<std::endl;
    
    //if (!read_sigmajets_from_file)
    //{
        //double tolerance=0.05, upper_tAA_0_limit=46.0, lower_tAA_0_limit = 30.0, upper_sumTpp_limit=0.61, lower_sumTpp_limit=0.03411;
      //  double tolerance=0.05, upper_tAA_0_limit=42.0, lower_tAA_0_limit = 33.0, upper_sumTpp_limit=0.5, lower_sumTpp_limit=0.03411;
      //  std::variant<InterpMultilinear<4, xsectval>, linear_interpolator, xsectval> 
      //      sigma_jets = calculate_spatial_sigma_jets_mf(tolerance, p_pdf, mand_s, kt02, jet_params, upper_tAA_0_limit, lower_tAA_0_limit, upper_sumTpp_limit, lower_sumTpp_limit);
        ////array<spatial,4> args{0.1, 0.3, 35.0, 40.0};
        ////std::cout<<sigma_jets.interp(args.begin())<<std::endl;
        ////std::cout<<pqcd::calculate_spatial_sigma_jet_mf(p_pdf, p_pdf, &mand_s, &kt02, &jet_params, &args[0], &args[1], &args[2], &args[3])<<std::endl;

        //double tolerance=0.25;//, upper_tAA_0_limit=46.0, lower_tAA_0_limit = 30.0, upper_sumTpp_limit=0.61, lower_sumTpp_limit=0.03411;
        //                                            1.35134e-186, 1.4788e-177, 1.02829e-215       4.24149, 3.609, 10.7237
        //const std::array<const double, 4>  lower_limits = {0.0, 0.0, 0.0, 30.0}, upper_limits = {4.5, 4.5, 11.0, 46.0};
        //const std::array<const double, 4>  lower_limits = {0.0, 0.0, 0.0, 28.4}, upper_limits = {3.6, 3.6, 8.8, 41.2};
        //const std::array<const double, 4>  lower_limits = {0.0, 0.0, 0.0, 30.5}, upper_limits = {3.6, 3.6, 8.8, 30.5};
        //double tolerance=0.05, upper_tAA_0_limit=42.0, lower_tAA_0_limit = 33.0, upper_sumTpp_limit=0.5, lower_sumTpp_limit=0.03411;
        //InterpMultilinear<5, xsectval> sigma_jets = calculate_spatial_sigma_jets_full(tolerance, p_pdf, p_pdf, mand_s, kt02, jet_params, lower_limits, upper_limits);
        //array<spatial,5> args{3.5, 4.3, 8.0, 39.0, 42.0}; //-304.782
        //std::cout<<sigma_jets.interp(args.begin())<<std::endl;

        //std::cout<<pqcd::calculate_spatial_sigma_jet_mf(p_pdf, p_pdf, &mand_s, &kt02, &jet_params, &args[0], &args[1], &args[2], &args[3])<<std::endl;
    //}    

    if (verbose) std::cout<<"Done!"<<std::endl;
    
    //std::ofstream log_file;
    //log_file.open("log2.dat", std::ios::out);



    auto cmpLambda = [](const io::Coll &lhs, const io::Coll &rhs) { return io::compET(lhs, rhs); };
    std::vector<std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)> > colls_scatterings(6, std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)>(cmpLambda));

    std::vector<histo_2d> jets(6, histo_2d({kt_bins, y_bins}));
    std::vector<histo_2d> dijets(6, histo_2d({kt_bins, y_bins}));
    std::mutex jets_mutex; 
    std::mutex dijets_mutex;

    std::vector<histo_1d> dETdy(6, histo_1d({y_bins}));
    std::vector<histo_1d> dEdy(6, histo_1d({y_bins}));
    std::vector<histo_1d> dNdy(6, histo_1d({y_bins}));
    std::vector<histo_1d> dNdET(6, histo_1d({et_bins}));
    std::mutex dETdy_mutex; 
    std::mutex dEdy_mutex;
    std::mutex dNdy_mutex;
    std::mutex dNdET_mutex;

    std::vector<histo_1d> dETdeta(6, histo_1d({y_bins}));
    std::vector<histo_1d> dEdeta(6, histo_1d({y_bins}));
    std::mutex dETdeta_mutex; 
    std::mutex dEdeta_mutex;

    std::array<std::ofstream, 6> total_energy;
    std::mutex total_energy_mutex; 
    auto max_energy = 0.5*(diff_params.A+diff_params.B)*sqrt_s;
    std::array<int64_t, 6> max_energy_broken{0,0,0,0,0,0};
 
    total_energy[0].open("total_energies_"+name_postfix+".dat");
    total_energy[1].open("total_energies_"+name_postfix+"_MC.dat");
    total_energy[2].open("total_energies_"+name_postfix+"_MC_ND.dat");
    total_energy[3].open("total_energies_"+name_postfix+"_SAT.dat");
    total_energy[4].open("total_energies_"+name_postfix+"_SAT_OL.dat");
    total_energy[5].open("total_energies_"+name_postfix+"_SAT_MC.dat");
    for (size_t i=0; i<6; i++)
    {
        total_energy[i] << "///Sum E_T Sum E" << std::endl;
    }
 
    std::vector<histo_1d> dETdb(6, histo_1d({b_bins}));
    std::vector<histo_1d> dEdb(6, histo_1d({b_bins}));
    std::mutex dETdb_mutex; 
    std::mutex dEdb_mutex;

    if (!read_events_from_file)
    {
        std::ofstream event_file;
        std::mutex event_file_mutex; 

        if (save_events)
        {
            event_file.open(event_file_name);

            event_file << "///AA b: ["<<b_min<<", "<<b_max<<"] fm"<<std::endl;
            event_file << "///"<<std::endl;
            event_file << "///Nucleus params:"<<std::endl;
            event_file << "///Nucleus r: ["<<rad_min<<", "<<rad_max<<"] fm"<<std::endl;
            event_file << "///Correct overlap bias: "<<nuc_params.correct_overlap_bias<<std::endl;
            event_file << "///Nucleon min_distance: "<<nuc_params.min_distance<<" fm"<<std::endl;
            event_file << "///Shift CMS: "<<nuc_params.shift_cms<<std::endl;
            event_file << "///"<<std::endl;
            event_file << "///Collision params:"<<std::endl;
            event_file << "///Proton width: "<<sqrt(proton_width_2)<<" fm"<<std::endl;
            event_file << "///Sqrt(s): "<<sqrt_s<<" GeV"<<std::endl;
            event_file << "///kt0: "<<kt0<<" GeV"<<std::endl;
            event_file << "///Projectile with npdfs: "<<diff_params.projectile_with_npdfs<<std::endl;
            event_file << "///Target with npdfs: "<<diff_params.target_with_npdfs<<std::endl;
            event_file << "///Spatial npdfs: "<<coll_params.spatial_pdfs<<std::endl;
            event_file << "///Projectile isoscalar: "<<diff_params.isoscalar_projectile<<std::endl;
            event_file << "///Target isoscalar: "<<diff_params.isoscalar_target<<std::endl;
            event_file << "///Normalization mode: "<<coll_params.normalize_to<<std::endl;
            event_file << "///Scale choice: "<<jet_params.scale_c<<std::endl;
            event_file << "///Scalar: "<<jet_params.scalar<<std::endl;
            event_file << "///"<<std::endl;
            event_file << "///End state params:"<<std::endl;
            event_file << "///Calculate end state: "<<coll_params.calculate_end_state<<std::endl;
            event_file << "///Power law: "<<power_law<<std::endl;
            event_file << "///"<<std::endl;
        }

        std::vector<uint64_t> event_indexes(desired_N_events);
        std::iota(event_indexes.begin(), event_indexes.end(), 0); //generates the list as {0,1,2,3,...}
        std::atomic<uint64_t> running_count{desired_N_events};

        std::signal(SIGINT, abort_handler);
        //std::signal(SIGTERM, term_handler);
        //std::signal(SIGSEGV, term_handler);

        std::set_terminate([](){
            std::cout << std::endl << "Unhandled exception" << std::endl;
            g_bug_bool = true;
        });

        try
        {
            std::find_if
            (
                std::execution::par, 
                event_indexes.begin(), 
                event_indexes.end(), 
                [&,&b_max=b_max,
                 &g_use_snpdfs=g_use_snpdfs,
                 &g_is_mom_cons_new=g_is_mom_cons_new,
                 &g_is_saturation=g_is_saturation,
                 &g_is_sat_overlap=g_is_sat_overlap,
                 &g_are_ns_depleted=g_are_ns_depleted,
                 &sigma_jet=sigma_jet,
                 &power_law=power_law,
                 &envelope_maximum=envelope_maximum,
                 &K_sat=K_sat,
                 &colls_scatterings=colls_scatterings
                ](const uint64_t index) 
                {
                    do //while (g_bug_bool)
                    {
                        uint32_t NColl = 0;
                        std::array<std::vector<nn_coll>, 6> binary_collisions;
                        std::array<std::vector<dijet_with_coords>, 6> filtered_scatterings;

                        //B^2 from a uniform distribution
                        spatial impact_parameter = 0;//sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min)); 

                        //std::shared_ptr<ars> radial_sampler;
                        //bool sampler_found = false;
                        //uint8_t s_index = index % samplersN;

                        //{
                        //    std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                        //    do //while sampler_found
                        //    {
                        //        if (!sampler_pointers[s_index]->locked())
                        //        {
                        //            radial_sampler = sampler_pointers[s_index];
                        //            sampler_found = true;
                        //        }
                        //        s_index = (s_index + 1) % samplersN;
                        //        //std::this_thread::sleep_for(std::chrono::milliseconds(10));
                        //    }
                        //    while (!sampler_found);
                        //    //std::cout<<std::endl<<s_index<<std::endl;
                        //}
                        //std::unique_lock<std::mutex> lock_rad(radial_sampler_mutex);
                        //auto radial_sampler = sampler_pointers[0];
                        std::vector<nucleon> pro, tar;

                        bool bugged, collided;
                        // "ball" diameter = distance at which two nucleons interact
                        const spatial d2 = 70.0/(M_PI*10.0); // in fm^2
                        do //while (!collided || bugged)
                        {
                            bugged = false;
                            collided = false;
                            try
                            {
                                auto [pro_dummy, tar_dummy] = calcs::generate_nuclei
                                (
                                    nuc_params, 
                                    sqrt_s, 
                                    impact_parameter, 
                                    eng, 
                                    radial_sampler, 
                                    read_nuclei_from_file, 
                                    verbose
                                );
                                pro = std::move(pro_dummy);
                                tar = std::move(tar_dummy);
                            }
                            catch(const std::exception& e)
                            {
                                std::cout << e.what() << " in main, trying again"<<std::endl;
                                bugged = true;
                            }
                            
                            for (auto A : pro)
                            {
                                for (auto B : tar)
                                {
                                    // "ball" diameter = distance at which two nucleons interact
                                    const spatial dij2 = A.calculate_bsquared(B);
                                    
                                    if (dij2 <= d2) //at least one collision
                                    {
                                        collided = true;
                                        continue;
                                    }
                                }
                                if (collided)
                                {
                                    continue;
                                }
                            }

                        } while (!collided || bugged);
                        
                        do //while (NColl<1)
                        {
                            binary_collisions[0].clear();
                            if (verbose) std::cout<<"impact_parameter: "<<impact_parameter<<std::endl;

                            calcs::collide_nuclei
                            (
                                g_use_snpdfs,
                                pro, 
                                tar, 
                                binary_collisions[0], 
                                sigma_jet,
                                unirand, 
                                eng, 
                                coll_params, 
                                jet_params,
                                kt0,
                                p_pdf,
                                power_law,
                                envelope_maximum,
                                verbose
                            );

                            /*if (average_spatial_taas)
                            {
                                //collide_nuclei_with_spatial_pdfs_averaging(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp);
                            }
                            else
                            {
                                //log_file << impact_parameter << std::endl;
                                //collide_nuclei_with_spatial_pdfs_factored(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp, proton_width_2, {impact_parameter,0,0}, log_file);
                                //collide_nuclei_with_spatial_pdfs_full(pro, tar, binary_collisions, sigma_jets, unirand, eng, coll_params, verbose, Tpp, proton_width_2, {impact_parameter,0,0}, log_file);
                                //log_file << std::endl;
                            }*/
                            NColl = static_cast<uint32_t>(binary_collisions[0].size());
                            if (NColl<1)
                            {
                                impact_parameter = 0;//sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min));
                                //const std::lock_guard<std::mutex> lock(radial_sampler_mutex);
                                std::tie(pro, tar) = calcs::generate_nuclei
                                (
                                    nuc_params, 
                                    sqrt_s, 
                                    impact_parameter, 
                                    eng, 
                                    radial_sampler, 
                                    read_nuclei_from_file, 
                                    verbose
                                );
                            }
                        } while (NColl<1);


                        for (size_t i=1; i<6; i++)
                        {
                            binary_collisions[i].assign(binary_collisions[0].begin(), binary_collisions[0].end());
                        }
                        std::array<momentum, 6> sum_ET{0,0,0,0,0,0};
                        std::array<momentum, 6> sum_ET_midrap{0,0,0,0,0,0};

                        if(!end_state_filtering && save_events)
                        {
                            const std::lock_guard<std::mutex> lock(event_file_mutex);
                            io::save_event(event_file, pro, tar, impact_parameter);
                        }
                        else if (end_state_filtering)
                        {
                            double max_overlap = K_sat;
                            std::array<momentum, 6> sum_E{0,0,0,0,0,0};
                            filter_end_state
                            (
                                binary_collisions[0], 
                                filtered_scatterings[0], 
                                unirand, 
                                eng, 
                                false, 
                                false,
                                false, 
                                false, 
                                false, 
                                sqrt_s,
                                max_overlap,
                                proton_width,
                                sigma_inel_for_glauber,
                                include_tata,
                                pro,
                                tar,
                                coll_params.Tpp
                            );
                            filter_end_state
                            (
                                binary_collisions[1], 
                                filtered_scatterings[1], 
                                unirand, 
                                eng, 
                                true, 
                                false,
                                false, 
                                true, 
                                false, 
                                sqrt_s, 
                                max_overlap,
                                proton_width,
                                sigma_inel_for_glauber,
                                include_tata,
                                pro,
                                tar,
                                coll_params.Tpp
                            );
                            filter_end_state
                            (
                                binary_collisions[2], 
                                filtered_scatterings[2], 
                                unirand, 
                                eng, 
                                true, 
                                false,
                                false, 
                                false, 
                                false, 
                                sqrt_s, 
                                max_overlap,
                                proton_width,
                                sigma_inel_for_glauber,
                                include_tata,
                                pro,
                                tar,
                                coll_params.Tpp
                            );
                            filter_end_state
                            (
                                binary_collisions[3], 
                                filtered_scatterings[3], 
                                unirand, 
                                eng, 
                                false, 
                                false,
                                true, 
                                false, 
                                false, 
                                sqrt_s, 
                                max_overlap,
                                proton_width,
                                sigma_inel_for_glauber,
                                include_tata,
                                pro,
                                tar,
                                coll_params.Tpp
                            );
                            filter_end_state
                            (
                                binary_collisions[4], 
                                filtered_scatterings[4], 
                                unirand, 
                                eng, 
                                false, 
                                false,
                                true, 
                                false, 
                                true, 
                                sqrt_s, 
                                max_overlap,
                                proton_width,
                                sigma_inel_for_glauber,
                                include_tata,
                                pro,
                                tar,
                                coll_params.Tpp
                            );
                            filter_end_state
                            (
                                binary_collisions[5], 
                                filtered_scatterings[5], 
                                unirand, 
                                eng, 
                                true,
                                false, 
                                true, 
                                true, 
                                false, 
                                sqrt_s, 
                                max_overlap,
                                proton_width,
                                sigma_inel_for_glauber,
                                include_tata,
                                pro,
                                tar,
                                coll_params.Tpp
                            );

                            std::array<std::vector<std::tuple<double, double> >, 6> new_jets;
                            std::array<std::vector<std::tuple<double, double> >, 6> new_dijets;
                            std::array<std::vector<std::tuple<double, double> >, 6> new_ET_y;
                            std::array<std::vector<std::tuple<double, double> >, 6> new_E_y;
                            std::array<std::vector<std::tuple<double, double> >, 6> new_N_y;
                            std::array<std::vector<std::tuple<double, double> >, 6> new_ET_eta;
                            std::array<std::vector<std::tuple<double, double> >, 6> new_E_eta;

                            for (size_t i=0; i<6; i++)
                            {
                                for (auto e_co : filtered_scatterings[i])
                                {
                                    auto e = e_co.dijet;

                                    new_jets[i].emplace_back(e.kt, e.y1);
                                    new_jets[i].emplace_back(e.kt, e.y2);
                                    
                                    new_dijets[i].emplace_back(e.kt, 0.5*(e.y1+e.y2));

                                    new_ET_y[i].emplace_back(e.y1, e.kt);
                                    new_ET_y[i].emplace_back(e.y2, e.kt);
                                    
                                    new_ET_eta[i].emplace_back(0.5*(e.y1+e.y2), 2*e.kt);
                                    
                                    new_E_y[i].emplace_back(e.y1, e.kt*cosh(e.y1));
                                    new_E_y[i].emplace_back(e.y2, e.kt*cosh(e.y2));
                                    
                                    new_N_y[i].emplace_back(e.y1, 1);
                                    new_N_y[i].emplace_back(e.y2, 1);
                                    
                                    new_E_eta[i].emplace_back(0.5*(e.y1+e.y2), e.kt*(cosh(e.y1) + cosh(e.y2)));
                                    
                                    sum_ET[i] += 2*e.kt;
                                    sum_E[i] += e.kt*(cosh(e.y1) + cosh(e.y2));

                                    if (e.y1 >= -0.5 && e.y1 <= 0.5)
                                    {
                                        sum_ET_midrap[i] += e.kt;
                                    }
                                    if (e.y2 >= -0.5 && e.y2 <= 0.5)
                                    {
                                        sum_ET_midrap[i] += e.kt;
                                    }
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(jets_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    jets[i].add(new_jets[i]);
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dijets_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dijets[i].add(new_dijets[i]);
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dETdy_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dETdy[i].add(new_ET_y[i]);
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dEdy_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dEdy[i].add(new_E_y[i]);
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dNdy_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dNdy[i].add(new_N_y[i]);
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dETdeta_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dETdeta[i].add(new_ET_eta[i]);
                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dEdeta_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dEdeta[i].add(new_E_eta[i]);
                                }
                            }
                            for (size_t i=0; i<6; i++)
                            {
                                {
                                    const std::lock_guard<std::mutex> lock(total_energy_mutex);
                                    total_energy[i] << sum_ET[i] << ' ' << sum_E[i] << std::endl;
                                }
                                
//                                if (sum_E[i] > max_energy)
//                                {
//                                    std::cout << std::endl 
//                                            << "Energy conservation violated! Total "<< ++max_energy_broken[i]
//                                            << " times this far in branch "<< i << std::endl
//                                            << "E_T = " << sum_ET[i] << ", E = " << sum_E[i]  << std::endl;
//                                }
                            }

                            {
                                const std::lock_guard<std::mutex> lock(dNdET_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dNdET[i].add(std::make_tuple(sum_ET[i], 1));
                                }
                            }
                            
                            {
                                const std::lock_guard<std::mutex> lock(dETdb_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dETdb[i].add(std::make_tuple(impact_parameter, sum_ET[i]));
                                }
                            }
                            {
                                const std::lock_guard<std::mutex> lock(dEdb_mutex);
                                for (size_t i=0; i<6; i++)
                                {
                                    dEdb[i].add(std::make_tuple(impact_parameter, sum_E[i]));
                                }
                            }

                            if (save_events)
                            {
                                const std::lock_guard<std::mutex> lock(event_file_mutex);
                                io::save_event(event_file, pro, tar, impact_parameter, filtered_scatterings[0]);
                            }

                            //filtered_scatterings.erase(filtered_scatterings.begin(), filtered_scatterings.end());
                        }
                        else
                        {
                            //binary_collisions.erase(binary_collisions.begin(), binary_collisions.end());
                        }

                        /*if (verbose || (nof_collisions%100)==0)
                        {
                            std::cout << std::endl << "A+A collided thus far: " << nof_collisions << ", of which events thus far: " << AA_events << std::endl << std::endl;
                        }*/

                        {
                            const std::lock_guard<std::mutex> lock(AA_events_mutex);
                            AA_events++;
                            if (AA_events % 100 == 0 )
                            {
                                std::cout <<"\rA+A collisions calculated: " << AA_events << std::flush;
                            }
                        }
                        //PrintThread{} <<"\rA+A collisions calculated: " << AA_events << std::flush;

                        uint32_t Npart=0;
                        for (auto &A : pro)
                        {
                            if (A.wounded)
                            {
                                Npart++;
                            }
                        }
                        for (auto &B : tar)
                        {
                            if (B.wounded)
                            {
                                Npart++;
                            }
                        }

                        {
                            const std::lock_guard<std::mutex> lock(colls_mutex);
                            for (size_t i=0; i<6; i++)
                            {
                                io::Coll coll(NColl, Npart, 2*filtered_scatterings[i].size(), impact_parameter, sum_ET[i]);
                                io::Coll coll_midrap(NColl, Npart, 2*filtered_scatterings[i].size(), impact_parameter, sum_ET_midrap[i]);
                                collisions_for_reporting[i].push_back(coll);
                                collisions_for_reporting_midrap[i].push_back(coll_midrap);

                                colls_scatterings[i].insert({coll, filtered_scatterings[i]});
                            }
                            
                            //io::save_single_coll(filtered_scatterings,unirand,eng,include_tata);
                        }
                        
                        //filtered_scatterings.clear();
                
                    } while (g_bug_bool);

                    bool ret_value = user_aborted;
                    return ret_value;
                }
            );
        }
        catch(const std::exception& e)
        {
            std::cout<<std::endl<<"Threw " << e.what() <<std::endl;
        }
        std::cout<<" ...done!"<<std::endl;
    }
    /*
    else
    {
        std::ifstream event_file;
        event_file.open(event_file_name);

        if (event_file.is_open())
        {
            std::string line;
            std::getline(event_file, line);
            //Skip comments at the start
            while (line[0] == '/')
            {
                std::getline(event_file, line);
            }

            do
            {
                std::getline(event_file, line);
                std::istringstream line_stream(line);
                spatial impact_parameter;
                line_stream >> impact_parameter;

                std::vector<nucleon> pro;
                std::getline(event_file, line);
                
                //Read projectile
                do 
                {
                    std::getline(event_file, line);
                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    coords co;
                    double x, y, z;

                    line_stream.ignore(256,'{');
                    line_stream >> x;
                    line_stream.ignore(256,',');
                    line_stream >> y;
                    line_stream.ignore(256,',');
                    line_stream >> z;
                    co = coords{x, y, z};

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    momentum mom;
                    line_stream >> mom;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool wounded;
                    line_stream >> wounded;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool is_neutron;
                    line_stream >> is_neutron;

                    pro.push_back(nucleon{co, mom, wounded, is_neutron});

                    std::getline(event_file, line);
                } while (line.back() != '}');
                pro.shrink_to_fit();

                std::getline(event_file, line);
                std::vector<nucleon> tar;
                std::getline(event_file, line);

                //Read target
                do
                {
                    std::getline(event_file, line);
                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    coords co;
                    double x, y, z;

                    line_stream.ignore(256,'{');
                    line_stream >> x;
                    line_stream.ignore(256,',');
                    line_stream >> y;
                    line_stream.ignore(256,',');
                    line_stream >> z;
                    co = coords{x, y, z};

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    momentum mom;
                    line_stream >> mom;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool wounded;
                    line_stream >> wounded;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    bool is_neutron;
                    line_stream >> is_neutron;

                    tar.push_back(nucleon{co, mom, wounded, is_neutron});

                    std::getline(event_file, line);
                } while (line.back() != '}');
                tar.shrink_to_fit();
                
                std::getline(event_file, line);
                std::getline(event_file, line);

                //Read filtered scatterings
                do
                {
                    std::getline(event_file, line);
                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    momentum kt;
                    line_stream >> kt;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    rapidity y1;
                    line_stream >> y1;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    rapidity y2;
                    line_stream >> y2;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id init1;
                    line_stream >> init1;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id init2;
                    line_stream >> init2;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id final1;
                    line_stream >> final1;

                    std::getline(event_file, line);
                    line_stream = std::istringstream(line);
                    particle_id final2;
                    line_stream >> final2;

                    filtered_scatterings.push_back(dijet_specs{kt, y1, y2, init1, init2, final1, final2});

                    std::getline(event_file, line);
                } while (line.back() != '}');
                filtered_scatterings.shrink_to_fit();

                std::getline(event_file, line);
                std::getline(event_file, line);

                nof_collisions++;
                AA_events++;
                momentum ET=0, E=0;
                std::vector<std::tuple<double, double> > new_jets;
                std::vector<std::tuple<double, double> > new_dijets;
                std::vector<std::tuple<double, double> > new_ET_y;
                std::vector<std::tuple<double, double> > new_E_y;
                std::vector<std::tuple<double, double> > new_ET_eta;
                std::vector<std::tuple<double, double> > new_E_eta;

                for (auto e : filtered_scatterings)
                {
                    new_jets.emplace_back(e.kt, e.y1);
                    new_jets.emplace_back(e.kt, e.y2);

                    new_dijets.emplace_back(e.kt, 0.5*(e.y1+e.y2));

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_ET_eta.emplace_back(0.5*(e.y1+e.y2), 2*e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));

                    new_E_eta.emplace_back(0.5*(e.y1+e.y2), e.kt*(cosh(e.y1) + cosh(e.y2)));

                    ET += 2*e.kt;
                    E += e.kt*(cosh(e.y1) + cosh(e.y2));
                }
 //               jets.add(new_jets);
 //               dijets.add(new_dijets);

                //dETdy.add(new_ET_y);
                //dEdy.add(new_E_y);

                //dETdeta.add(new_ET_eta);
                //dEdeta.add(new_E_eta);

 //               dETdb.add(std::make_tuple(impact_parameter, ET));
 //               dEdb.add(std::make_tuple(impact_parameter, E));

                filtered_scatterings.erase(filtered_scatterings.begin(), filtered_scatterings.end());

                std::cout<<"\rA+A collisions read thus far: " << nof_collisions << std::flush;

                std::getline(event_file, line);
            } while (line[0] == '{');
        }
        else
        {
            std::cout<<"ERROR READING EVENTS"<<std::endl;
        }
    }*/


    io::print_histos
    (
        name_postfix,
        jets,
        dijets,
        dETdy,
        dEdy,
        dNdy,
        dNdET,
        dETdeta,
        dEdeta,
        dETdb,
        dEdb,
        dijet_norm,
        AA_events
    );

    //log_file.close();

    uint nBins = 18;
    double binsLow[] = {0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 0.0, 0.0, 0.0, 0.0};
    double binsHigh[] = {0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 0.05, 0.1, 0.8, 1.0};
    std::ofstream glauber_report_file;
    std::array<std::string, 6> namesss{"g_report_"+name_postfix+".dat",
                                       "g_report_"+name_postfix+"_MC.dat",
                                       "g_report_"+name_postfix+"_MC_ND.dat",
                                       "g_report_"+name_postfix+"_SAT.dat",
                                       "g_report_"+name_postfix+"_SAT_OL.dat",
                                       "g_report_"+name_postfix+"_SAT_MC.dat"};
    std::array<std::string, 6> namesss_midrap{"g_report_midrap_"+name_postfix+".dat",
                                              "g_report_midrap_"+name_postfix+"_MC.dat",
                                              "g_report_midrap_"+name_postfix+"_MC_ND.dat",
                                              "g_report_midrap_"+name_postfix+"_SAT.dat",
                                              "g_report_midrap_"+name_postfix+"_SAT_OL.dat",
                                              "g_report_midrap_"+name_postfix+"_SAT_MC.dat"};
    for (size_t i=0; i<6; i++)
    {
        glauber_report_file.open(namesss[i], std::ios::out);
        io::mc_glauber_style_report(collisions_for_reporting[i], sigma_inel_for_glauber, desired_N_events, nBins, binsLow, binsHigh, glauber_report_file);
        glauber_report_file.close();
        glauber_report_file.open(namesss_midrap[i], std::ios::out);
        io::mc_glauber_style_report(collisions_for_reporting_midrap[i], sigma_inel_for_glauber, desired_N_events, nBins, binsLow, binsHigh, glauber_report_file);
        glauber_report_file.close();
    }

    
    double centLow = 0.0;
    double centHigh = 0.05;
    std::array<std::string, 6> names_pfs{name_postfix+".dat",
                                         name_postfix+"_MC.dat",
                                         name_postfix+"_MC_ND.dat",
                                         name_postfix+"_SAT.dat",
                                         name_postfix+"_SAT_OL.dat",
                                         name_postfix+"_SAT_MC.dat"};
    for (size_t i=0; i<6; i++)
    {
        histo_1d dETdy_by_cent{y_bins};
        histo_1d dEdy_by_cent{y_bins};
        std::vector<std::tuple<double, double> > new_ET_y;
        std::vector<std::tuple<double, double> > new_E_y;

        size_t N_evts_tot = colls_scatterings[i].size();
        // Make sure that no rounding downwards.
        double eps = 0.1/N_evts_tot;

        size_t lower_ind = static_cast<size_t>(centLow*N_evts_tot+eps);
        size_t upper_ind = static_cast<size_t>(centHigh*N_evts_tot+eps);

        for (auto it = colls_scatterings[i].crbegin(); lower_ind<upper_ind; ++it, lower_ind++)
        {
            for (auto e_co : it->second)
            {
                auto e = e_co.dijet;

                new_ET_y.emplace_back(e.y1, e.kt);
                new_ET_y.emplace_back(e.y2, e.kt);
                
                new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
            }
        }

        dETdy_by_cent.add(new_ET_y);
        dEdy_by_cent.add(new_E_y);

        io::print_1d_histo
        (
            dETdy_by_cent,
            "dETdy_0_5_"+names_pfs[i], 
            1.0/ (N_evts_tot*(centHigh-centLow)),
            false
        );
        io::print_1d_histo
        (
            dEdy_by_cent, 
            "dEdy_0_5_"+names_pfs[i], 
            1.0/ (N_evts_tot*(centHigh-centLow)),
            false
        );
    }

    return 0;
}