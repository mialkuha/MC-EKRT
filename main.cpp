//Copyright (c) 2022 Mikko Kuha 

#include <algorithm>
#include <atomic>
#include <execution>
#include <csignal>
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

auto find_sigma_jet_cutoff
(
    momentum &kt02, 
    const momentum &mand_s, 
    const xsectval &target, 
    std::shared_ptr<LHAPDF::GridPDF> p_p_pdf, 
    const pqcd::sigma_jet_params &jet_params, 
    const bool &verbose=true
) noexcept -> void
{
    xsectval sigma_jet=0.0;
    kt02 = 2.0;

    auto difference_to_target = [&](const momentum &_kt02)
    {
        return pqcd::calculate_sigma_jet(p_p_pdf, &mand_s, &_kt02, jet_params) - target;
    };

    helpers::secant_method(&kt02, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<kt02<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return;
}
auto check_and_place_circle_among_others
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

struct dijet_with_ns
{
    dijet_specs dijet;
    nucleon * pro_nucleon;
    nucleon * tar_nucleon;
    dijet_with_ns(dijet_specs dijet_, nucleon * pro_nucleon_, nucleon * tar_nucleon_)
        : dijet(std::move(dijet_)), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_) { }
};

auto throw_location_for_dijet //TODO
(
    const dijet_with_ns &cand
) noexcept -> std::tuple<double, double>
{
    return std::make_tuple(0, 0);
}

//Checks whether the saturation criterium allows the candidate to be added TODO
auto filter_collisions_saturation
(
    std::vector<dijet_with_ns> &candidates, 
    std::vector<dijet_specs> &final_candidates,
    const double & maximum_overlap,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> void
{
    std::vector<std::tuple<double, double, double, uint16_t> > final_circles;
    final_circles.reserve(candidates.size());
    final_candidates.reserve(final_candidates.size()+candidates.size());

    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the biggest kt is first
        [](dijet_with_ns &s1, dijet_with_ns &s2) { return (s1.dijet.kt > s2.dijet.kt); });

    for (auto & cand : candidates)
    {
        auto [cand_x, cand_y] = throw_location_for_dijet(cand);
        auto cand_circle = std::make_tuple(cand_x, cand_y, 1/cand.dijet.kt, 1);
        if (check_and_place_circle_among_others(std::move(cand_circle), final_circles, maximum_overlap, unirand, random_generator))
        {
            final_candidates.emplace_back(std::move(cand.dijet));
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
    colls_with_ns(momentum kt_, rapidity y1_, rapidity y2_, nucleon * pro_nucleon_, nucleon * tar_nucleon_, dijet_specs * dijet_)
        : kt(kt_), y1(y1_), y2(y2_), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_), dijet(dijet_) { }
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

    binary_collisions.clear();
    //for (auto &col : binary_collisions)
    //{
    //    col.dijets.erase(col.dijets.begin(), col.dijets.end());
    //}
}

auto filter_end_state
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_specs> &filtered_scatterings,
    std::uniform_real_distribution<double> &unirand,
    std::shared_ptr<std::mt19937> random_generator,
    const bool mom_cons = false,
    const bool saturation = false,
    const bool deplete_nucleons = false,
    const momentum sqrt_s = 0,
    const double & maximum_overlap = 2.0
) noexcept -> void
{
    std::vector<dijet_with_ns> candidates;
    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess
    
    if (mom_cons) 
    {
        filter_collisions_MC(binary_collisions, candidates, sqrt_s, deplete_nucleons);
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
        filter_collisions_saturation(candidates, filtered_scatterings, maximum_overlap, unirand, random_generator);
    }
    else
    {
        for (auto &candidate : candidates) //all of the candidates are passed to filtered
        {
            filtered_scatterings.push_back(candidate.dijet);
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
    bool g_is_pa;
    bool g_is_aa;
    bool g_use_npdfs;
    bool g_use_snpdfs;
    bool g_is_mom_cons;
    bool g_is_mom_cons_new;
    bool g_are_ns_depleted;
    bool g_is_saturation;
    std::string name_postfix;

    switch (argv[1][0])
    {
    case '0':
        name_postfix = "_pA_1m_mb_PDF";
        g_is_pa = false;
        g_is_aa = true;
        g_use_npdfs = false;
        g_use_snpdfs = false;
        g_is_mom_cons = false;
        g_is_mom_cons_new = false;
        g_are_ns_depleted = false;
        g_is_saturation = false;
        break;
    case '1':
        name_postfix = "_pA_1m_mb_nPDF";
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
        name_postfix = "_pA_1m_mb_snPDF";
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
        name_postfix = "_pA_1m_mb_PDF_MC";
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
        name_postfix = "_pA_1m_mb_nPDF_MC";
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
        name_postfix = "_pA_1m_mb_snPDF_MC";
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
        name_postfix = "_pA_1m_mb_PDF_MC_ND";
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
        name_postfix = "_pA_1m_mb_nPDF_MC_ND";
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
        name_postfix = "_pA_1m_mb_snPDF_MC_ND";
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
    }

    //A lot of printing
    bool verbose = false;
    
    //General parameters for the simulation
    const bool    read_nuclei_from_file    = false, 
                  read_events_from_file    = false,
                  read_sigmajets_from_file = true,
                  end_state_filtering      = true, 
                  save_events              = false/*, 
                  average_spatial_taas     = false*/;
    std::string   event_file_name = "event_log"+name_postfix+".dat";
    uint32_t      desired_N_events      = 1000000000,
                  AA_events             = 0;
    const spatial b_min                 = 0,
                  b_max                 = 20;
    auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand{0.0, 1.0};

    const spatial rad_min=0, rad_max=30;
    const std::function<double(const double&)> rad_pdf{[](const double & x)
    {
        return x*x/(1+exp((x-6.624)/0.549));
    }};
    auto radial_sampler{std::make_shared<ars>(rad_pdf, rad_min, rad_max)};
    do //while (radial_sampler->is_adaptive())
    {
        radial_sampler->throw_one(*eng);
    } while (radial_sampler->is_adaptive());

    nucleus_generator::nucleus_params nuc_params = 
    {
        /* .NA=                   */(g_is_aa)? 208u : 1u, //Pb 
        /* .ZA=                   */(g_is_aa)? 82u : 1u, //Pb 
        /* .NB=                   */208, 
        /* .ZB=                   */82, 
        /* .min_distance=         */0.4, 
        /* .shift_cms=            */true, 
        /* .correct_overlap_bias= */true
    };
    const momentum sqrt_s                 = 5020;//GeV

    if (!read_events_from_file)
    {

        std::vector<uint64_t> event_indexes(desired_N_events);
        std::iota(event_indexes.begin(), event_indexes.end(), 0); //generates the list as {0,1,2,3,...}

        std::signal(SIGINT, abort_handler);
        std::signal(SIGTERM, term_handler);
        std::signal(SIGSEGV, term_handler);

        std::set_terminate([](){
            std::cout << std::endl << "Unhandled exception" << std::endl;
            g_bug_bool = true;
        });

        std::mutex event_file_mutex; 
        try
        {
            std::find_if
            (
                std::execution::par, 
                event_indexes.begin(), 
                event_indexes.end(), 
                [&](const uint64_t index) 
                {
                bug_bagk:
                    spatial impact_parameter = sqrt(b_min*b_min + unirand(*eng)*(b_max*b_max-b_min*b_min)); 

                    auto [pro, tar] = calcs::generate_nuclei
                    (
                        nuc_params, 
                        sqrt_s, 
                        impact_parameter, 
                        eng, 
                        radial_sampler, 
                        read_nuclei_from_file, 
                        verbose
                    );
                    {
                        const std::lock_guard<std::mutex> lock(event_file_mutex);

                        if (g_bug_bool)
                        {
                            goto bug_bagk;
                        }

                        AA_events++;
                        if (AA_events % 1000 == 0 )
                        {
                            std::cout <<"\rA+A coords calculated: " << AA_events << std::flush;
                        }
                    }

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

    return 0;
}
