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
    // PUBLIC ATTRIBUTES //////////////////////////////////////////////////////////////////////////////////////////////
    // technical flags
    bool calculate_end_state{true};        // false=no jets will be calculated, only collisions
    bool calculate_tata{true};             // if TA(x)TA(x) should be calculated for each dijet (takes time)
    bool deplete_nucleons{false};          // probably always false here, kept for testing purposes for now
    bool end_state_filtering{true};        // false=no any filtering in the end state
    bool is_AA{true};                      // if the run is nucleus-nucleus
    bool is_pA{false};                     // if the run is proton-nucleus
    bool is_pp{false};                     // if the run is proton-proton
    bool nPDFs{true};                      // if there should be nuclear corrections in nucleus' PDFs
    bool MC_Glauber{true};                 // if true, nucleon-nucleon collisions are decided like with hard spheres
    bool mom_cons{true};                   // if the momentum should be conserved
    bool read_sigmajets_from_file{false};  // if the sigma_jets are precalculated
    bool reduce_nucleon_energies{false};   // momentum conservation by reducing nucleon energies after each collision. Breaks factorization
    bool proton_width_static{false};       // if the proton width should be a static value or depend on sqrt(s)
    bool saturation{true};                 // if EKRT saturation is enforced
    bool save_endstate_jets{true};         // if all the jets should be saved in binary (see jet_reader.cpp)
    bool save_events_plaintext{false};     // if all the jets should be saved in plaintext
    bool sigma_inel_from_sigma_jet{true};  // is sigma_inel static or from eikonal model
    bool snPDFs{false};                    // if the used nPDF should be spatially dependent or average
    bool verbose{false};                   // lots of printing all around
    // simulation parameters 
    std::string name{"example_name"};      // name of the run, affects output filenames and such 
    uint_fast32_t desired_N_events{500};   // how many events should be simulated
    double b_max{20.0};                    // (fm) maximum of the impact parameter
    double b_min{0.0};                     // (fm) minimum of the impact parameter 
    double K_factor{1.0};                  // pQCD K-factor to account for the higher order corrections
    double kt0{1.0};                       // (GeV) jet p_T lower cutoff
    double kt02{1.0};                      // (GeV^2) p_T lower cutoff squared
    double M_factor{2.5};                  // saturation criterion for jets i, j: d_ij < (1/p_Ti + 1/p_Tj)/M_sat
    double nn_min_dist{0.4};               // (fm) minimum distance between the nucleons in a nucleus
    double mand_s{5020*5020};              // (GeV^2) mandelstam s for the hard process 
    double proton_width{0.573};            // (fm) proton width
    double proton_width_2{0.573*0.573};    // (fm^2) proton width squared
    double rad_max{30.0};                  // (fm) nucleons' maximum distance from the center of the nucleus
    double rad_min{0.0};                   // (fm) nucleons' minimum distance from the center of the nucleus
    double sigma_inel{70.0};               // (mb) inelastic cross section of the event
    double sqrt_s{5020};                   // (GeV) sqrt(s) for the hard process 


    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
    /**
     * @brief Construct a new mcaa object, some of the parameters
     * are taken from a parameter file. See "params_template" for 
     * an example file.
     * 
     * @param initfile input parameter file name
     */
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
     * @brief Runs the simulator with the given parameters. The reporting is
     * also built in this function, so that if one would like to change the 
     * outputs for example, this is the function to modify. This is done to
     * keep in memory only the things that are actually needed.
     * 
     */
    auto run() -> void;


// PRIVATE ////////////////////////////////////////////////////////////////////////////////////////////////////////////
private:
    // PRIVATE ATTRIBUTES /////////////////////////////////////////////////////////////////////////////////////////////
    std::shared_ptr<LHAPDF::GridPDF> p_pdf{};     // pointer to the LHAPDF PDF-object of the proton PDF
    std::function<double(const double&)> Tpp{nullptr}; // nucleon overlap function
    std::function<double(const double&)> rad_pdf{nullptr}; // radial nucleon density function of the nucleus
    double power_law{2.0};                        // parameter for the envelope function: sigma_jet < A*pT^(-power_law)
    pqcd::diff_sigma::params diff_params;         // the struct of pQCD event parameters
    pqcd::sigma_jet_params jet_params;            // the struct of sigma_jet parameters
    nucleus_generator::nucleus_params nuc_params; // the struct of parameters for the nucleus generation


    // PRIVATE STRUCTS ////////////////////////////////////////////////////////////////////////////////////////////////
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

    
    // PRIVATE METHODS ////////////////////////////////////////////////////////////////////////////////////////////////
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
    static auto check_and_place_circle_among_others
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
     * @param pro Projectile nucleus' nucleon information
     * @param tar Target nucleus' nucleon information
     */
    auto filter_end_state
    (
    /*uint_fast64_t &nall,
    uint_fast64_t &nmc,
    uint_fast64_t &nsat,*/
        std::vector<nn_coll> &binary_collisions, 
        std::vector<dijet_with_coords> &filtered_scatterings,
        std::shared_ptr<std::mt19937> random_generator,
        const std::vector<nucleon> &pro = {}, 
        const std::vector<nucleon> &tar = {}
    ) noexcept -> void;
};

#endif // MCAA_HPP
