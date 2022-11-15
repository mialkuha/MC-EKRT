//Copyright (c) 2022 Mikko Kuha

#ifndef MCAA_HPP
#define MCAA_HPP

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

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "LHAPDF/GridPDF.h"
#pragma GCC diagnostic pop

#include "ars.hpp"
#include "generic_helpers.hpp"
#include "high_level_calcs.hpp"
#include "histo.hpp"
#include "io_helpers.hpp"
#include "nucleus_generator.hpp"
#include "pqcd.hpp"
#include "typedefs.hpp"

class mcaa
{
public:

    mcaa(const std::string &initfile);

    /**
     * @brief Find k_T^2 cutoff that fits sigma_jet to a given 
     * target number using secant method.
     * 
     * @param pt02 k_T^2 cutoff will be saved to this (GeV^2)
     * @param target the target value for sigma_jet (mb)
     * @param verbose if true, prints all the intermediate results
     * @return double the found k_T^2 cutoff (GeV^2)
     */
    auto fit_sigma_jet_pt0_cutoff
    (
        double &pt02, 
        const double &target,
        const bool &verbose=true
    ) noexcept -> double;

    /**
     * @brief Find scale choice Q/k_T that fits sigma_jet to a given 
     * target number, for a given target k_T^2 cutoff, using secant method.
     * 
     * @param pt02 The target k_T^2 cutoff value (GeV^2)
     * @param target the target value for sigma_jet (mb)
     * @param verbose if true, prints all the intermediate results
     * @return double the found scale choice Q/k_T
     */
    auto find_sigma_jet_cutoff_Q
    (
        const double &pt02, 
        const double &target,
        const bool &verbose=true
    ) noexcept -> double;

    /**
     * @brief TODOTODO
     * 
     */
    auto run() -> void;

// PRIVATE ///////////////////////////////////////////////////////////
private:
    // FIELDS /////////////////////////////
    bool verbose{false};                      //lots of printing all around
    bool mom_cons{false};                     //if the momentum should be conserved
    bool saturation{false};                   //if EKRT saturation is enforced
    bool deplete_nucleons{false};             //TODO explanation 
    bool is_AA{false};             //TODO explanation 
    bool is_pA{false};             //TODO explanation 
    bool is_pp{false};             //TODO explanation 
    bool nPDFs{false};             //TODO explanation 
    bool snPDFs{false};             //TODO explanation 
    bool MC_Glauber{false};             //TODO explanation 
    std::string name{""};             //TODO explanation 
    uint_fast32_t desired_N_events{1};
    double K_factor{1.0};                     //factor for the saturation criterion
    double M_factor{2.0};                     //factor for the saturation criterion
    double kt0{1.0};                          //k_T lower cutoff for end state jets
    double kt02{1.0};                         //k_T lower cutoff squared
    double sqrt_s{5020};                      //sqrt(s) for the hard process (GeV^2)
    double mand_s{5020*5020};                 //mandelstam s for the hard process (GeV^2)
    double proton_width{1};               //proton width (fm)
    double proton_width_2{1};               //proton width (fm)
    double b_min{0};           //TODO explanation 
    double b_max{0};           //TODO explanation 
    std::shared_ptr<LHAPDF::GridPDF> p_pdf{}; //pointer to the LHAPDF PDF-object of the proton PDF
    pqcd::diff_sigma::params diff_params;
    pqcd::sigma_jet_params jet_params; //the struct of jet parameters
    std::function<double(const double&)> Tpp{nullptr}; //nucleon overlap function

    // STRUCTS ////////////////////////////
    /**
     * @brief struct for the end state filtering functions. In the constructor 
     * calculates the formation time t0 = y/max(u,t) = (1/Q)*(E/Q) (in the collider frame) 
     * for the jets.
     * 
     */
    struct dijet_with_ns
    {
        dijet_specs dijet;
        nucleon * pro_nucleon;
        nucleon * tar_nucleon;
        double t01;
        double t02;
        double t0; //Larger of the formation times := the formation time of the process
        dijet_with_ns(dijet_specs dijet_, nucleon * pro_nucleon_, nucleon * tar_nucleon_)
            : dijet(std::move(dijet_)), pro_nucleon(pro_nucleon_), tar_nucleon(tar_nucleon_)
        {
            auto y1 = dijet.y1;
            auto y2 = dijet.y2;
            auto max_t_u = dijet.kt*(1+std::exp(-std::abs(y1-y2)));
            t01 = std::cosh(std::abs(y1))/max_t_u;
            t02 = std::cosh(std::abs(y2))/max_t_u;
            t0 = std::max(t01, t02);
        }
    };

    
    // FUNCTIONS /////////////////////////////
    /**
    * @brief Checks if the candidate circle overlaps with the previous circles. If
    * it does not, add it to the collection of accepted circles.
    * 
    * @param cand_circle (x,y,r,n), where x and y are coordinates (fm), r is 
    * radius (fm) and n is for supporting, say, two circles overlapping but no 
    * more (not implemented).
    * @param final_circles A collection of the previously accepted circles.
    * @return true The circle was accepted and added to the collection.
    * @return false The circle was not accepted.
    */
    auto check_and_place_circle_among_others
    (
        std::tuple<double, double, double, uint_fast16_t> cand_circle, 
        std::vector<std::tuple<double, double, double, uint_fast16_t> > &final_circles
    ) noexcept -> bool;

    /**
     * @brief Sample the coordinates for a given dijet from the Gaussian distribution
     * that is the product of the two Gaussian distributions that are associated with
     * the mother nucleons.
     * 
     * @param dijet The dijet in question
     * @param normal_dist Normal distribution object
     * @param random_generator Random generator object
     * @return std::tuple<double, double, double> (x,y,z=0) coordinates of the dijet (fm)
     */
    auto throw_location_for_dijet
    (
        const dijet_with_ns &dijet,
        std::normal_distribution<double> &normal_dist, 
        std::shared_ptr<std::mt19937> random_generator
    ) noexcept -> std::tuple<double, double, double>;

    /**
     * @brief Takes the unphysical final state of the hard process and filters it with the 
     * momentum conservation criterion and EKRT saturation criterion.
     * 
     * @param binary_collisions The set of (unphysical) pQCD final states
     * @param filtered_scatterings The output, fully filtered dijets
     * @param random_generator Random generator object
     * @param calculate_tata If TA(x)TA(x) should be calculated for each dijet (takes time)
     * @param pro Projectile nucleus' nucleon information
     * @param tar Target nucleus' nucleon information
     */
    auto filter_end_state
    (
        std::vector<nn_coll> &binary_collisions, 
        std::vector<dijet_with_coords> &filtered_scatterings,
        std::shared_ptr<std::mt19937> random_generator,
        const bool calculate_tata = true,
        const std::vector<nucleon> &pro = {}, 
        const std::vector<nucleon> &tar = {}
    ) noexcept -> void;
};

#endif // MCAA_HPP
