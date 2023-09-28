//Copyright (c) 2023 Mikko Kuha

#include "mcaa.hpp"

mcaa::mcaa
(
    const std::string &initfile
) 
: diff_params(false,false,false,false,false,false,1,0,1,1u,1u,1u,1u),
  jet_params(this->diff_params,pqcd::scaled_from_kt,1,false),
  nuc_params(1u,1u,1u,1u,1.0,false,false)
{
    auto 
    [
        p_name,
        p_sjet_name,
        p_cent_name,
        p_n_events,
        p_b_max,
        p_b_min,
        p_sqrt_s,
        p_K_factor,
        p_kt0,
        p_proton_width,
        p_sigma_inel,
        p_sigma_inel_AA,
        p_T_AA_0_for_snpdfs,
        p_spatial_cutoff,
        p_A,
        p_ZA,
        p_M_factor,
        p_nn_min_dist,
        p_nuclear_RA,
        p_nuclear_d,
        p_rad_max,
        p_rad_min,
        p_is_AA,
        p_is_pA,
        p_is_pp,
        p_is_mc_glauber,
        p_read_sigmajets_from_file,
        p_proton_width_static,
        p_sigma_inel_from_sigma_jet,
        p_AA_inel_same_as_NN,
        p_only_protons,
        p_use_npdfs,
        p_use_snpdfs,
        p_snpdfs_linear,
        p_calculate_spatial_cutoff,
        p_calculate_end_state,
        p_calculate_tata,
        p_save_endstate_jets,
        p_end_state_filtering,
        p_is_mom_cons,
        p_is_saturation,
        p_is_sat_y_dep,
        p_is_pt_ordering,
        p_is_t03_ordering
    ] = io::read_conf(initfile);

    this->calculate_end_state        = p_calculate_end_state; 
    this->calculate_tata             = p_calculate_tata; 
    this->end_state_filtering        = p_end_state_filtering; 
    this->is_AA                      = p_is_AA; 
    this->is_pA                      = p_is_pA; 
    this->is_pp                      = p_is_pp; 
    this->nPDFs                      = p_use_npdfs; 
    this->MC_Glauber                 = p_is_mc_glauber; 
    this->mom_cons                   = p_is_mom_cons; 
    this->read_sigmajets_from_file   = p_read_sigmajets_from_file; 
    this->proton_width_static        = p_proton_width_static; 
    this->saturation                 = p_is_saturation; 
    this->save_endstate_jets         = p_save_endstate_jets;
    this->sigma_inel_from_sigma_jet  = p_sigma_inel_from_sigma_jet; 
    this->AA_inel_same_as_NN         = p_AA_inel_same_as_NN; 
    this->only_protons               = p_only_protons; 
    this->snPDFs                     = p_use_snpdfs;
    this->snPDFs_linear              = p_snpdfs_linear;
    this->is_sat_y_dep               = p_is_sat_y_dep;
    this->pt_ordering                = p_is_pt_ordering;
    this->t03_ordering               = p_is_t03_ordering;
    this->name                       = p_name;
    this->sigmajet_filename          = p_sjet_name;
    this->centrality_filename        = p_cent_name;
    this->desired_N_events           = p_n_events;
    this->A                          = p_A;
    this->ZA                         = p_ZA;
    this->b_max                      = p_b_max;
    this->b_min                      = p_b_min;
    this->K_factor                   = p_K_factor;
    this->kt0                        = p_kt0;
    this->M_factor                   = p_M_factor;
    this->nn_min_dist                = p_nn_min_dist;
    this->nuclear_RA                 = p_nuclear_RA;
    this->nuclear_d                  = p_nuclear_d;
    this->rad_max                    = p_rad_max;
    this->rad_min                    = p_rad_min;
    this->sigma_inel                 = p_sigma_inel;
    this->sigma_inel_AA              = p_sigma_inel_AA;
    this->sqrt_s                     = p_sqrt_s;
    this->T_AA_0                     = p_T_AA_0_for_snpdfs;
    this->spatial_cutoff             = p_spatial_cutoff;

    this->verbose                    = false;
    this->kt02                       = p_kt0*p_kt0;
    this->mand_s                     = p_sqrt_s*p_sqrt_s;


    if (!this->AA_inel_same_as_NN && almostEqual(p_sigma_inel_AA, 0.0)) // triggering sigma from parametrization
    {
        auto s = this->mand_s;
        // fit by COMPETE https://doi.org/10.1103/PhysRevLett.89.201801
        auto sigmatot = 42.6*std::pow(s,-0.46) -33.4*std::pow(s,-0.545) +0.307*std::pow(std::log(s/29.1),2) +35.5;
        // fit by TOTEM https://doi.org/10.1140/epjc/s10052-019-6567-0
        auto sigmael = -1.617*std::log(s) +0.1359*std::pow(std::log(s),2) +11.84;

        this->sigma_inel_AA = sigmatot - sigmael;
        std::cout<<"Triggering sigma_inel from parametrization: "<<this->sigma_inel_AA<<std::endl;
    }

    if (this->proton_width_static)
    {
        this->proton_width           = p_proton_width;
        this->proton_width_2         = p_proton_width*p_proton_width;
    }
    else
    {
        this->proton_width_2 = (4.9 + 4*0.06*std::log(this->sqrt_s/90.0))/(FMGEV*FMGEV); //C. A. Flett, PhD thesis, U. Liverpool, 2021
        this->proton_width = std::sqrt(this->proton_width_2);
    }
    std::cout<<"Proton width: "<<this->proton_width<<" fm"<<std::endl;

    this->p_pdf = std::make_shared<LHAPDF::GridPDF>("CT14lo", 0);

    this->Tpp =  std::function<double(const double&)>({[proton_width_2=this->proton_width_2](const double &bsquared)
    {
        return exp(-bsquared / (4 * proton_width_2)) / (40 * M_PI * proton_width_2); // 1/fm² = mb/fm² * 1/mb = 0.1 * 1/mb
    }});

    if (p_use_snpdfs && p_snpdfs_linear && p_calculate_spatial_cutoff)
    {
        auto A = this->A;
        auto RA = this->nuclear_RA;
        auto d = this->nuclear_d;
        auto tp0 = 1.0 / (2.0 * M_PI * this->proton_width_2);
        auto n0 = 0.75 * (A/(M_PI*std::pow(RA,3.0))) / (1.0 + M_PI*M_PI*std::pow(d/RA,2.0));
        auto ta0 = 2.0 * n0 * d * std::log(1.0 + std::exp(RA/d));
        this->spatial_cutoff = tp0 / (3.0 * ta0); 
        std::cout<<"Calculated spatial cutoff: (1+cT) >= "<<this->spatial_cutoff<<std::endl;
    }

    this->rad_pdf = std::function<double(const double&)>({[RA=this->nuclear_RA, d=this->nuclear_d](const double & x)
    {
        return x*x/(1+exp((x-RA)/d));
    }});

    // Parameter for the envelope function:
    // sigma_jet < A*pT^(-power_law)
    // Larger number == faster simulation, but the
    // inequality must hold for the results to be correct.
    // These are just some values, experiment freely,
    // the simulation will tell if the inequality breaks.
    if (!this->nPDFs) //free proton PDFs 
    {
        this->power_law = 2.8;
    }
    else if (!this->snPDFs) //nPDFs
    {
        this->power_law = 2.2;
    }
    else //spatial nPDFs
    {
        this->power_law = 2.0;
    }
    
    this->diff_params = pqcd::diff_sigma::params(
    /*projectile_with_npdfs=    */(p_is_AA && p_use_npdfs),
    /*target_with_npdfs=        */(!p_is_pp && p_use_npdfs),
    /*isoscalar_projectile=     */false,
    /*isoscalar_target=         */false,
    /*npdfs_spatial=            */p_use_snpdfs,
    /*only_protons=             */p_only_protons,
    /*npdf_setnumber=           */1,
    /*spatial_cutoff=           */this->spatial_cutoff,
    /*K_factor=                 */K_factor,
    /*A=                        */(p_is_AA)? p_A : 1u,  //Pb 
    /*B=                        */(p_is_pp)? 1u : p_A,  //Pb
    /*ZA=                       */(p_is_AA)? p_ZA : 1u, //Pb
    /*ZB=                       */(p_is_pp)? 1u : p_ZA  //Pb
    /*p_n_pdf=                  */
    /*rA_spatial=               */
    /*rB_spatial=               */);

    this->jet_params = pqcd::sigma_jet_params(
    /*d_params=                 */this->diff_params,
    /*scale_choice=             */pqcd::scaled_from_kt,
    /*scalar=                   */1.0,
    /*use_ses=                  */false);

    this->nuc_params = nucleus_generator::nucleus_params(
    /*A=                        */(p_is_AA)? p_A : 1u,  //Pb 
    /*ZA=                       */(p_is_AA)? p_ZA : 1u, //Pb
    /*B=                        */(p_is_pp)? 1u : p_A,  //Pb
    /*ZB=                       */(p_is_pp)? 1u : p_ZA, //Pb
    /*min_distance=             */p_nn_min_dist, 
    /*shift_cms=                */true, 
    /*correct_overlap_bias=     */true);
}

auto mcaa::fit_sigma_jet_pt0_cutoff
(
    double &pt02, 
    const double &target, 
    const bool &verbose
) noexcept -> double
{
    double sigma_jet=0.0;
    pt02 = 2.0;

    auto difference_to_target = [&](const double &_kt02)
    {
        return pqcd::calculate_sigma_jet(this->p_pdf, &(this->mand_s), &_kt02, this->jet_params) - target;
    };

    helpers::secant_method(&pt02, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<sqrt(pt02)<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return pt02;
}

auto mcaa::find_sigma_jet_cutoff_Q
(
    const double &pt02, 
    const double &target,
    const bool &verbose
) noexcept -> double
{
    double sigma_jet=0.0;
    double scalar = 1.0;

    auto difference_to_target = [&](const double &scalar_)
    {
        double kt02dummy = 4.0;
        auto jet_params_ = pqcd::sigma_jet_params(
        /*d_params=                 */this->jet_params.d_params,
        /*scale_choice=             */this->jet_params.scale_c,
        /*scalar=                   */scalar_,
        /*use_ses=                  */this->jet_params.use_ses);
        auto kt02_ = this->fit_sigma_jet_pt0_cutoff(kt02dummy, target);
        return kt02_ - pt02;
    };

    helpers::secant_method(&scalar, difference_to_target, 1e-3, &sigma_jet);

    if (verbose) std::cout<<scalar<<' '<<sigma_jet+target<<' '<<target<<std::endl;
    
    return scalar;
}

auto mcaa::check_and_place_circle_among_others
(
    std::tuple<double, double, double, uint_fast16_t, double, double> cand_circle, 
    std::vector<std::tuple<double, double, double, uint_fast16_t, double, double> > &final_circles,
    uint_fast16_t is_sat_y_dep
) noexcept -> bool
{
    auto & [cand_x, cand_y, cand_r, cand_overlap, cand_y1, cand_y2] = cand_circle;
    
    for (auto & [circ_x, circ_y, circ_r, circ_overlap, y1, y2] : final_circles)
    {
        //If the any of the circles radii overlap with the candidate, just return false
        if ( pow((circ_x-cand_x),2) + pow((circ_y-cand_y),2) < pow((circ_r+cand_r),2) )
        {
            if (is_sat_y_dep == 1) // Henry's old idea
            {
                auto miy = std::min(y1,y2);
                auto may = std::max(y1,y2);
                auto mciy = std::min(cand_y1,cand_y2);
                auto mcay = std::max(cand_y1,cand_y2);
                
                if (!(miy>mcay || may<mciy))
                {
                    return false;
                }
            }
            else if (is_sat_y_dep == 2) // all jets have cones
            {
                double delta = 0.5;
                if ((std::abs(y1-cand_y1) <= delta)
                 || (std::abs(y1-cand_y2) <= delta)
                 || (std::abs(y2-cand_y1) <= delta)
                 || (std::abs(y2-cand_y2) <= delta))
                {
                    return false;
                }
            }
            else if (is_sat_y_dep == 3) // both dijets have cones
            {
                double delta = 0.5;
                auto dy = 0.5*(y1+y2);
                auto cand_dy = 0.5*(cand_y1+cand_y2);

                if ((std::abs(dy-cand_dy) <= delta))
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }

    //None overlapped, emplace the circle among others and return true
    final_circles.emplace_back(std::move(cand_circle));
        
    return true;
}

auto mcaa::throw_location_for_dijet
(
    const dijet_with_ns &dijet,
    std::normal_distribution<double> &normal_dist, 
    std::shared_ptr<std::mt19937> random_generator
) noexcept -> std::tuple<double, double, double>
{
    auto param = std::normal_distribution<double>::param_type{0., this->proton_width};
    auto dx = normal_dist(*random_generator,param);
    auto dy = normal_dist(*random_generator,param);
    auto dz = 0.0;

    auto retx = 0.5*(dijet.pro_nucleon->co.x + dijet.tar_nucleon->co.x + M_SQRT2*dx);
    auto rety = 0.5*(dijet.pro_nucleon->co.y + dijet.tar_nucleon->co.y + M_SQRT2*dy);

    return std::make_tuple(retx, rety, dz);
}

auto mcaa::sqrtalphas(double Q) noexcept -> double
{
    auto alphas = this->p_pdf->alphasQ(Q);
    return std::sqrt(alphas);
}

auto mcaa::filter_end_state
(
    std::vector<nn_coll> &binary_collisions, 
    std::vector<dijet_with_coords> &filtered_scatterings,
    std::vector<dijet_with_coords> &filtered_scatterings_nomc,
    std::vector<dijet_with_coords> &filtered_scatterings_nosat,
    std::vector<dijet_with_coords> &filtered_scatterings_neither, 
    std::shared_ptr<std::mt19937> random_generator,
    const std::vector<nucleon> &pro, 
    const std::vector<nucleon> &tar
) noexcept -> void
{
    std::vector<dijet_with_ns> candidates;
    candidates.reserve(binary_collisions.size()*10); //10 events on average is just an overhead guess

    std::unordered_map<nucleon*, double> x1s;
    std::unordered_map<nucleon*, double> x2s;
    

    std::normal_distribution<double> normal_dist(0,0);

    for (auto &col : binary_collisions)
    {
        for (auto &dij : col.dijets)
        {
            candidates.emplace_back(std::move(dij), col.projectile, col.target, this->pt_ordering, this->t03_ordering);
        }
    }
    candidates.shrink_to_fit();
    
    std::vector<std::tuple<double, double, double, uint_fast16_t, double, double> > final_circles;
    final_circles.reserve(candidates.size());
    filtered_scatterings.reserve(filtered_scatterings.size()+candidates.size());
    
    std::sort(candidates.begin(), candidates.end(), //Sort the candidate events so that the one with the smallest t0 is first
              [](dijet_with_ns &s1, dijet_with_ns &s2) { return (s1.t0 < s2.t0); });
    
    double kt;
    double y1, y2;
    double i_x1_sum_to_be, i_x2_sum_to_be;
    double tata = 0.0;


    for (auto & cand : candidates)
    {
        bool mc_broke = false;
        bool sat_broke = false;

        kt = cand.dijet.kt;
        y1 = cand.dijet.y1;
        y2 = cand.dijet.y2;

        if (this->mom_cons)
        {
            //MOMENTUM CONSERVATION

            auto x1 = (kt / this->sqrt_s) * (exp(y1) + exp(y2));
            auto x2 = (kt / this->sqrt_s) * (exp(-y1) + exp(-y2));

            i_x1_sum_to_be = x1;
            i_x2_sum_to_be = x2;

            auto i_x1_sum = x1s.find(cand.pro_nucleon);
            if (i_x1_sum != x1s.end())
            {
                i_x1_sum_to_be = i_x1_sum->second + x1;

                if (i_x1_sum_to_be > 1.0) //Energy budget broken --> discard
                {
                    mc_broke = true;
                }
            }

            auto i_x2_sum = x2s.find(cand.tar_nucleon);
            if (i_x2_sum != x2s.end())
            {
                i_x2_sum_to_be = i_x2_sum->second + x2;

                if (i_x2_sum_to_be > 1.0) //Energy budget broken --> discard
                {
                    mc_broke = true;
                }
            }
        }
        
        auto [cand_x, cand_y, cand_z] = this->throw_location_for_dijet(cand, normal_dist, random_generator);

        if (this->saturation)
        {
            //SATURATION

            double upstairs = 1.0;
            if (this->is_sat_y_dep == 4)
            {
                upstairs = this->sqrtalphas(cand.dijet.kt);
            }
            else if  (this->is_sat_y_dep == 5)
            {
                double dy = y1-y2;
                double dsigma = std::pow((1 + 2.0*std::cosh(dy))/(2 + 2.0*std::cosh(dy)), 1.5);
                upstairs = this->sqrtalphas(cand.dijet.kt)*dsigma;
            }

            auto cand_circle = std::make_tuple(cand_x, cand_y, upstairs/(cand.dijet.kt*this->M_factor*FMGEV), 1, y1, y2);

            if (!mcaa::check_and_place_circle_among_others(std::move(cand_circle), final_circles, this->is_sat_y_dep))
            {
                sat_broke = true;
            }
        }

        if (this->calculate_tata)
        {
            nucleon dummy{coords{cand_x, cand_y, cand_z}, 0};
            tata = calcs::calculate_sum_tpp(dummy, pro, this->Tpp) * calcs::calculate_sum_tpp(dummy, tar, this->Tpp);
        }

        if (this->mom_cons && !mc_broke)
        {
            x1s.insert_or_assign(cand.pro_nucleon, i_x1_sum_to_be);
            x2s.insert_or_assign(cand.tar_nucleon, i_x2_sum_to_be);
        }

        if (mc_broke && sat_broke)
        {
            filtered_scatterings_neither.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
        }
        else if (mc_broke)
        {
            filtered_scatterings_nomc.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
            filtered_scatterings_neither.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
        }
        else if (sat_broke)
        {
            filtered_scatterings_nosat.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
            filtered_scatterings_neither.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
        }
        else
        {
            filtered_scatterings.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
            filtered_scatterings_nomc.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
            filtered_scatterings_nosat.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
            filtered_scatterings_neither.push_back({cand.dijet, coords{cand_x, cand_y, cand_z}, cand.t01, cand.t02, tata});
        }
    }

    binary_collisions.clear();
}

//--------------------------------------------//
// THESE ARE A HACK TO SAVE THE WORK THIS FAR //
// IN THE CASE OF SUDDEN TERMINATION          //
//--------------------------------------------//
volatile std::atomic_bool user_aborted = false;
volatile std::atomic_bool g_bug_bool = false;
void abort_handler(int num)
{
    user_aborted = true;
    std::cout<<std::endl<<"Abort called "<<num<<std::endl;
}

auto mcaa::run() -> void
{
    std::cout<<std::endl<<"Doing the run "<<this->name<<std::endl;

    if (this->verbose) std::cout<<"Initializing..."<<std::flush;
    
    uint_fast32_t AA_events_done = 0;
    std::mutex AA_events_mutex;

    //auto eng = std::make_shared<std::mt19937>(static_cast<ulong>(1));
    auto eng_shared = std::make_shared<std::mt19937>(static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_real_distribution<double> unirand{0.0, 1.0};

    auto radial_sampler{std::make_shared<ars>(rad_pdf, rad_min, rad_max)};
    //The adaptive algorithm in the sampler is not thread-safe,
    //so to run the program multithreaded let's first saturate the sampler
    do //while (radial_sampler->is_adaptive()) 
    {
        radial_sampler->throw_one(*eng_shared);
    } while (radial_sampler->is_adaptive());

    std::vector<io::Coll> collisions_for_reporting;/////
    std::vector<io::Coll> collisions_for_reporting_midrap;/////
    std::vector<io::Coll> collisions_for_reporting_nomc;/////
    std::vector<io::Coll> collisions_for_reporting_midrap_nomc;/////
    std::vector<io::Coll> collisions_for_reporting_nosat;/////
    std::vector<io::Coll> collisions_for_reporting_midrap_nosat;/////
    std::vector<io::Coll> collisions_for_reporting_neither;/////
    std::vector<io::Coll> collisions_for_reporting_midrap_neither;/////
    std::mutex colls_mutex; 
    std::vector<double> kt_bins{helpers::loglinspace(kt0, 100, 16u)};
    auto kt_v_dummy = helpers::loglinspace(200, sqrt_s/2.0, 5u);
    for (auto ktdummy : kt_v_dummy)
    {
        kt_bins.push_back(ktdummy);
    }
    
    double ylim = static_cast<double>(log(this->sqrt_s / this->kt0));
    const std::vector<double> y_bins{helpers::linspace(-ylim, ylim, 40u)};
    const std::vector<double> b_bins{helpers::linspace(this->b_min, this->b_max, 21u)};
    std::vector<double> et_bins{helpers::loglinspace(2*kt0, 30000, 31u)};
    
    std::vector<double> Rs_as_vector, cs_as_vector;

    //if linear snPDFs
    if (this->snPDFs && (this->T_AA_0 == 0.0) && this->snPDFs_linear)
    {
        std::cout<<"Calculating T_AA(0)"<<std::endl;

        this->T_AA_0 = calcs::calculate_T_AA_0
        (
            this->nuc_params,
            1e-5,
            this->Tpp,
            radial_sampler, 
            verbose
        );

        std::cout<<"Calculated T_AA(0) = "<< this->T_AA_0 <<std::endl;
    } //if exponential snPDFs
    else if (this->snPDFs && (this->T_AA_0 == 0.0))
    {
        std::cout<<"Calculating R_A - c_A table"<<std::endl;

        auto [R_A_table, c_A_table] = calcs::calculate_R_c_table
        (
            this->nuc_params,
            1e-5,
            this->Tpp,
            radial_sampler, 
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

        auto c_A_from_R = linear_interpolator(Rs_as_vector, cs_as_vector);
    
        this->jet_params = pqcd::sigma_jet_params(
        /*d_params=                 */this->jet_params.d_params,
        /*scale_choice=             */this->jet_params.scale_c,
        /*scalar=                   */this->jet_params.scalar,
        /*use_ses=                  */this->jet_params.use_ses,
        /*snPDFs_linear=            */this->snPDFs_linear,
        /*c_A_func=                 */c_A_from_R);
    }

    auto
    [
        dijet_norm,
        sigma_jet,
        mand_s_array,
        env_func,
        sqrt_s_array
    ] = 
    calcs::prepare_sigma_jets
    (
        this->read_sigmajets_from_file,
        this->p_pdf, 
        this->mand_s,
        this->sqrt_s,
        this->kt02, 
        this->kt0,
        this->power_law,
        this->jet_params,
        this->T_AA_0,
        this->sigmajet_filename
    );

    if (this->sigma_inel_from_sigma_jet && !this->snPDFs)
    {

        auto dummy = std::get<double>(sigma_jet) / (4 * M_PI * this->proton_width_2);
        dummy = dummy / 10; //mb -> fm²

        this->sigma_inel = (4 * M_PI * this->proton_width_2) * (M_EULER + std::log(dummy) + gsl_sf_expint_E1(dummy)) * 10; //1 mb = 0.1 fm^2
        std::cout<<"sigma_inel = "<<this->sigma_inel<<std::endl;
    }
    else if (this->sigma_inel_from_sigma_jet)
    {
        auto dummy = dijet_norm / (4 * M_PI * this->proton_width_2);
        dummy = dummy / 10; //mb -> fm²

        this->sigma_inel = (4 * M_PI * this->proton_width_2) * (M_EULER + std::log(dummy) + gsl_sf_expint_E1(dummy)) * 10; //1 mb = 0.1 fm^2
        std::cout<<"sigma_inel = "<<this->sigma_inel<<std::endl;
    }

    if (this->AA_inel_same_as_NN)
    {
        this->sigma_inel_AA = this->sigma_inel;
    }

    bool use_nn_b2_max = true; //TODO GET RID OF MAGIC NUMBERS
    double tpp_min = 1e-8; //TODO GET RID OF MAGIC NUMBERS

    double nn_b2_max = 1.0;
    if (use_nn_b2_max)
    {
        double dummy = 0; //not needed here
        auto difference_to_target = [tpp_min, tpp=(this->Tpp)](const double &b2)
        {
            return tpp_min - tpp(b2);
        };
        helpers::secant_method(&nn_b2_max, difference_to_target, tpp_min, &dummy);
        if (this->verbose)
        {
            std::cout<<"Found b2max = "<<nn_b2_max<<", Tpp(b2max) = "<<this->Tpp(nn_b2_max)<<std::endl;
        }
    }

    if (this->verbose) std::cout<<"Done!"<<std::endl;

    AA_collision_params coll_params
    {
    /*mc_glauber_mode=          */this->MC_Glauber,
    /*pp_scattering=            */this->is_pp,
    /*pA_scattering=            */this->is_pA,
    /*spatial_pdfs=             */this->snPDFs,
    /*calculate_end_state=      */this->calculate_end_state,
    /*use_nn_b2_max=            */use_nn_b2_max,
    /*sigma_inel=               */this->sigma_inel,
    /*Tpp=                      */this->Tpp,
    /*normalize_to=             */B2_normalization_mode::inelastic,
    /*sqrt_s=                   */this->sqrt_s,
    /*energy_threshold=         */this->kt0,
    /*nn_b2_max=                */nn_b2_max,
    /*T_AA_0=                   */this->T_AA_0
    };

    auto cmpLambda = [](const io::Coll &lhs, const io::Coll &rhs) { return io::compET(lhs, rhs); };
    std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)> colls_scatterings(cmpLambda);
    std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)> colls_scatterings_nomc(cmpLambda);
    std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)> colls_scatterings_nosat(cmpLambda);
    std::map<io::Coll, std::vector<dijet_with_coords>, decltype(cmpLambda)> colls_scatterings_neither(cmpLambda);

    // auto cmpLambda = [](const std::tuple<double, double> &lhs, const std::tuple<double, double> &rhs) { return std::get<0>(lhs) < std::get<0>(rhs); }; /////
    // std::map<io::Coll, std::tuple<double, double>, decltype(cmpLambda)> colls_sumet_midet(cmpLambda); /////

//    OBSERVABLES TO BE SAVED
//
//    histo_2d jets{kt_bins, y_bins};
//    histo_2d dijets{kt_bins, y_bins};
//    std::mutex jets_mutex; 
//    std::mutex dijets_mutex;
//
//    histo_1d dETdy{y_bins};
//    histo_1d dEdy{y_bins};
//    histo_1d dNdy{y_bins};
//    histo_1d dNdET{et_bins};
//    std::mutex dETdy_mutex; 
//    std::mutex dEdy_mutex;
//    std::mutex dNdy_mutex;
//    std::mutex dNdET_mutex;
//
//    histo_1d dETdeta{y_bins};
//    histo_1d dEdeta{y_bins};
//    std::mutex dETdeta_mutex; 
//    std::mutex dEdeta_mutex;
// 
//    histo_1d dETdb{b_bins};
//    histo_1d dEdb{b_bins};
//    std::mutex dETdb_mutex; 
//    std::mutex dEdb_mutex;


    std::ofstream total_energy;/////
    std::mutex total_energy_mutex;/////
 
    total_energy.open("total_energies_"+this->name+".dat");/////
    total_energy << "///Sum E_T Sum E" << std::endl;/////

    std::vector<uint_fast64_t> event_indexes(this->desired_N_events);
    std::iota(event_indexes.begin(), event_indexes.end(), 0); //generates the list as {0,1,2,3,...}
    std::atomic<uint_fast64_t> running_count{desired_N_events};


    //--------------------------------------------//
    // THESE ARE A HACK TO SAVE THE WORK THIS FAR //
    // IN THE CASE OF A UNEXPECTED TERMINATION    //
    //--------------------------------------------//
    std::signal(SIGINT, abort_handler);
    std::set_terminate([](){
        std::cout << std::endl << "Unhandled exception" << std::endl;
        //g_bug_bool = true;
    });

    try
    {
        //----------------------------------------------------------------//
        // THIS IS THE PARALLELLIZED MAIN LOOP FOR ALL OF THE EVENTS      //
        // INCLUDES NUCLEUS GENERATION, COLLISION CALCULATIONS, FILTERING //
        // AND OBSERVABLE CALCULATIONS                                    //
        //----------------------------------------------------------------//
        
        #pragma omp parallel
        {
            auto eng = std::make_shared<std::mt19937>(static_cast<ulong>((omp_get_thread_num() + 1))*static_cast<ulong>(std::chrono::system_clock::now().time_since_epoch().count()));

            #pragma omp for
            for (auto it = event_indexes.begin(); it < event_indexes.end(); it++) 
            {
                if (user_aborted)
                {
                    continue;
                }
                do //while (g_bug_bool)
                {
                    uint_fast32_t NColl = 0;
                    std::vector<nn_coll> binary_collisions;
                    std::vector<dijet_with_coords> filtered_scatterings;
                    std::vector<dijet_with_coords> filtered_scatterings_nomc;
                    std::vector<dijet_with_coords> filtered_scatterings_nosat;
                    std::vector<dijet_with_coords> filtered_scatterings_neither;
                    double impact_parameter;
                    double b_min2 = this->b_min*this->b_min;
                    double b_max2 = this->b_max*this->b_max;

                    std::vector<nucleon> pro, tar;
                    uint_fast16_t times_discarded = 0;
                    
                    //Keep generating nuclei until the triggering condition is met
                    bool bugged, triggered;
                    // "ball" diameter = distance at which two nucleons interact in MC Glauber
                    const double d2 = this->sigma_inel_AA/(M_PI*10.0); // in fm^2
                    do //while (!triggered || bugged)  
                    {
                        //B^2 from a uniform distribution
                        impact_parameter = sqrt(b_min2 + unirand(*eng)*(b_max2-b_min2));
                
                        times_discarded++;
                        if (times_discarded > 1000)
                        {
                            std::cout<<std::endl<<"Generated nuclei discarded over 1000 times. "
                                        <<"Check impact parameters and/or collsion probabilities."
                                        <<std::endl;
                            times_discarded = 0;
                        }

                        bugged = false;
                        triggered = false;
                        try
                        {
                            auto [pro_dummy, tar_dummy] = calcs::generate_nuclei
                            (
                                this->nuc_params, 
                                this->sqrt_s, 
                                impact_parameter, 
                                eng, 
                                radial_sampler, 
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
                                const double dij2 = A.calculate_bsquared(B);

                                if (dij2 <= d2) // triggering condition
                                {
                                    triggered = true;
                                    continue;
                                }
                            }
                            if (triggered)
                            {
                                continue;
                            }
                        }

                    } while (!triggered || bugged);

                    //Demand at least one hard scattering
                    do //while (NColl<1)
                    {
                        binary_collisions.clear();
                        if (verbose) std::cout<<"impact_parameter: "<<impact_parameter<<std::endl;

                        calcs::collide_nuclei
                        (
                            pro, 
                            tar, 
                            binary_collisions, 
                            sigma_jet,
                            unirand, 
                            eng, 
                            coll_params, 
                            this->jet_params,
                            this->kt0,
                            this->p_pdf,
                            this->power_law,
                            env_func,
                            verbose,
                            this->M_factor,
                            this->proton_width,
                            this->is_sat_y_dep
                        );
                        
                        NColl = static_cast<uint_fast32_t>(binary_collisions.size());
                    } while (NColl<1);
                    

                    double sum_ET = 0;
                    double sum_ET_midrap = 0;
                    double sum_ET_nomc = 0;
                    double sum_ET_midrap_nomc = 0;
                    double sum_ET_nosat = 0;
                    double sum_ET_midrap_nosat = 0;
                    double sum_ET_neither = 0;
                    double sum_ET_midrap_neither = 0;

                    if (end_state_filtering)
                    {
                        double sum_E = 0;
                        this->filter_end_state
                        (
                            binary_collisions, 
                            filtered_scatterings,
                            filtered_scatterings_nomc,
                            filtered_scatterings_nosat,
                            filtered_scatterings_neither, 
                            eng,
                            pro,
                            tar
                        );

                        // std::vector<std::tuple<double, double> > new_jets;
                        // std::vector<std::tuple<double, double> > new_dijets;
                        // std::vector<std::tuple<double, double> > new_ET_y;
                        // std::vector<std::tuple<double, double> > new_E_y;
                        // std::vector<std::tuple<double, double> > new_N_y;
                        // std::vector<std::tuple<double, double> > new_ET_eta;
                        // std::vector<std::tuple<double, double> > new_E_eta;
                        
                        for (auto e_co : filtered_scatterings)
                        {
                            auto e = e_co.dijet;

                            // new_jets.emplace_back(e.kt, e.y1);
                            // new_jets.emplace_back(e.kt, e.y2);
                            
                            // new_dijets.emplace_back(e.kt, 0.5*(e.y1+e.y2));
                            
                            // new_ET_y.emplace_back(e.y1, e.kt);
                            // new_ET_y.emplace_back(e.y2, e.kt);
                            
                            // new_ET_eta.emplace_back(0.5*(e.y1+e.y2), 2*e.kt);
                            
                            // new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                            // new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                            
                            // new_N_y.emplace_back(e.y1, 1);
                            // new_N_y.emplace_back(e.y2, 1);
                            
                            // new_E_eta.emplace_back(0.5*(e.y1+e.y2), e.kt*(cosh(e.y1) + cosh(e.y2)));
                            
                            sum_ET += 2*e.kt;
                            sum_E += e.kt*(cosh(e.y1) + cosh(e.y2));/////

                            if (e.y1 >= -0.5 && e.y1 <= 0.5)
                            {
                                sum_ET_midrap += e.kt;
                            }
                            if (e.y2 >= -0.5 && e.y2 <= 0.5)
                            {
                                sum_ET_midrap += e.kt;
                            }
                        }
                        for (auto e_co : filtered_scatterings_nomc)
                        {
                            auto e = e_co.dijet;
                            
                            sum_ET_nomc += 2*e.kt;

                            if (e.y1 >= -0.5 && e.y1 <= 0.5)
                            {
                                sum_ET_midrap_nomc += e.kt;
                            }
                            if (e.y2 >= -0.5 && e.y2 <= 0.5)
                            {
                                sum_ET_midrap_nomc += e.kt;
                            }
                        }
                        for (auto e_co : filtered_scatterings_nosat)
                        {
                            auto e = e_co.dijet;
                            
                            sum_ET_nosat += 2*e.kt;

                            if (e.y1 >= -0.5 && e.y1 <= 0.5)
                            {
                                sum_ET_midrap_nosat += e.kt;
                            }
                            if (e.y2 >= -0.5 && e.y2 <= 0.5)
                            {
                                sum_ET_midrap_nosat += e.kt;
                            }
                        }
                        for (auto e_co : filtered_scatterings_neither)
                        {
                            auto e = e_co.dijet;
                            
                            sum_ET_neither += 2*e.kt;

                            if (e.y1 >= -0.5 && e.y1 <= 0.5)
                            {
                                sum_ET_midrap_neither += e.kt;
                            }
                            if (e.y2 >= -0.5 && e.y2 <= 0.5)
                            {
                                sum_ET_midrap_neither += e.kt;
                            }
                        }

                        // {
                        //    const std::lock_guard<std::mutex> lock(jets_mutex);
                        //    jets.add(new_jets);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dijets_mutex);
                        //    dijets.add(new_dijets);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dETdy_mutex);
                        //    dETdy.add(new_ET_y);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dEdy_mutex);
                        //    dEdy.add(new_E_y);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dNdy_mutex);
                        //    dNdy.add(new_N_y);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dETdeta_mutex);
                        //    dETdeta.add(new_ET_eta);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dEdeta_mutex);
                        //    dEdeta.add(new_E_eta);
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dNdET_mutex);
                        //    dNdET.add(std::make_tuple(sum_ET, 1));
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dETdb_mutex);
                        //    dETdb.add(std::make_tuple(impact_parameter, sum_ET));
                        // }
                        
                        // {
                        //    const std::lock_guard<std::mutex> lock(dEdb_mutex);
                        //    dEdb.add(std::make_tuple(impact_parameter, sum_E));
                        // }
                        
                        {/////
                            const std::lock_guard<std::mutex> lock(total_energy_mutex);/////
                            total_energy << sum_ET << ' ' << sum_E << std::endl;/////
                        }/////
                    }

                    {
                        const std::lock_guard<std::mutex> lock(AA_events_mutex);
                        AA_events_done++;
                        if (AA_events_done % 100 == 0 )
                        {
                            std::cout <<"\rA+A collisions calculated: " << AA_events_done << std::flush;
                        }
                    }

                    uint_fast32_t Npart=0;
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
                        io::Coll coll(NColl, Npart, 2*filtered_scatterings.size(), impact_parameter, sum_ET);
                        io::Coll coll_midrap(NColl, Npart, 2*filtered_scatterings.size(), impact_parameter, sum_ET_midrap);
                        io::Coll coll_nomc(NColl, Npart, 2*filtered_scatterings_nomc.size(), impact_parameter, sum_ET_nomc);
                        io::Coll coll_midrap_nomc(NColl, Npart, 2*filtered_scatterings_nomc.size(), impact_parameter, sum_ET_midrap_nomc);
                        io::Coll coll_nosat(NColl, Npart, 2*filtered_scatterings_nosat.size(), impact_parameter, sum_ET_nosat);
                        io::Coll coll_midrap_nosat(NColl, Npart, 2*filtered_scatterings_nosat.size(), impact_parameter, sum_ET_midrap_nosat);
                        io::Coll coll_neither(NColl, Npart, 2*filtered_scatterings_neither.size(), impact_parameter, sum_ET_neither);
                        io::Coll coll_midrap_neither(NColl, Npart, 2*filtered_scatterings_neither.size(), impact_parameter, sum_ET_midrap_neither);
                        collisions_for_reporting.push_back(coll);
                        collisions_for_reporting_midrap.push_back(coll_midrap);
                        collisions_for_reporting_nomc.push_back(coll_nomc);
                        collisions_for_reporting_midrap_nomc.push_back(coll_midrap_nomc);
                        collisions_for_reporting_nosat.push_back(coll_nosat);
                        collisions_for_reporting_midrap_nosat.push_back(coll_midrap_nosat);
                        collisions_for_reporting_neither.push_back(coll_neither);
                        collisions_for_reporting_midrap_neither.push_back(coll_midrap_neither);
                        
                        if (save_endstate_jets)
                        {
                            colls_scatterings.insert({coll, filtered_scatterings});
                            colls_scatterings_nomc.insert({coll_nomc, filtered_scatterings_nomc});
                            colls_scatterings_nosat.insert({coll_nosat, filtered_scatterings_nosat});
                            colls_scatterings_neither.insert({coll_neither, filtered_scatterings_neither});
                        }
                    }
                } while (g_bug_bool);
            }
        }
    }
    catch(const std::exception& e)
    {
        std::cout<<std::endl<<"Threw " << e.what() <<std::endl;
    }
    std::cout<<" ...done!"<<std::endl;

    

    //io::print_histos
    //(
    //    name_postfix,
    //    jets,
    //    dijets,
    //    dETdy,
    //    dEdy,
    //    dNdy,
    //    dNdET,
    //    dETdeta,
    //    dEdeta,
    //    dETdb,
    //    dEdb,
    //    dijet_norm,
    //    AA_events_done
    //);
    
    uint_fast8_t nBins = 18;
    double binsLow[] = {0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 0.0, 0.0, 0.0, 0.0};
    double binsHigh[] = {0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 0.05, 0.1, 0.8, 1.0};
    std::ofstream glauber_report_file;
    std::string g_name{"g_report_"+this->name+".dat"};
    std::string g_name_midrap{"g_report_midrap_"+this->name+".dat"};
    std::string g_name_nomc{"g_report_"+this->name+"_nomc.dat"};
    std::string g_name_midrap_nomc{"g_report_midrap_"+this->name+"_nomc.dat"};
    std::string g_name_nosat{"g_report_"+this->name+"_nosat.dat"};
    std::string g_name_midrap_nosat{"g_report_midrap_"+this->name+"_nosat.dat"};
    std::string g_name_neither{"g_report_"+this->name+"_neither.dat"};
    std::string g_name_midrap_neither{"g_report_midrap_"+this->name+"_neither.dat"};
                                              
    glauber_report_file.open(g_name, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting, this->sigma_inel, collisions_for_reporting.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_midrap, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_midrap, this->sigma_inel, collisions_for_reporting_midrap.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_nomc, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_nomc, this->sigma_inel, collisions_for_reporting_nomc.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_midrap_nomc, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_midrap_nomc, this->sigma_inel, collisions_for_reporting_midrap_nomc.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_nosat, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_nosat, this->sigma_inel, collisions_for_reporting_nosat.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_midrap_nosat, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_midrap_nosat, this->sigma_inel, collisions_for_reporting_midrap_nosat.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_neither, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_neither, this->sigma_inel, collisions_for_reporting_neither.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();
    glauber_report_file.open(g_name_midrap_neither, std::ios::out);
    io::mc_glauber_style_report(collisions_for_reporting_midrap_neither, this->sigma_inel, collisions_for_reporting_midrap_neither.size(), nBins, binsLow, binsHigh, glauber_report_file);
    glauber_report_file.close();

    if (save_endstate_jets)
    {
        std::vector<std::tuple<double,double> > centBins;
        std::ifstream cents_input(this->centrality_filename);
        if (!cents_input.is_open())
        {
            std::cout<<"Could not open "<<this->centrality_filename<<std::endl;
            centBins.push_back(std::tuple<double,double>{0.0, 1.0});
            return;
        }
        std::string line;
        for (; std::getline(cents_input, line);)
        {
            std::stringstream ss(line);
            std::string token;
            std::getline(ss, token, ',');
            double centLow = std::stod(token);
            std::getline(ss, token, ',');
            double centHigh = std::stod(token);
            centBins.push_back(std::tuple<double,double>{centLow, centHigh});
        }
        cents_input.close();

        std::string name_pfs{this->name+".dat"};

        std::ofstream jet_file;

        for (auto [centLow, centHigh] : centBins)
        {
            std::stringstream jetsname{""};
            jetsname<<"jets_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            jet_file.open(jetsname.str(), std::ios::out | std::ios::binary);

            histo_1d dETdy_by_cent{y_bins};
            histo_1d dEdy_by_cent{y_bins};
            std::vector<std::tuple<double, double> > new_ET_y;
            std::vector<std::tuple<double, double> > new_E_y;

            uint_fast64_t N_evts_tot = colls_scatterings.size();
            // Make sure that no rounding downwards.
            double eps = 0.1/static_cast<double>(N_evts_tot);

            uint_fast64_t lower_ind = static_cast<uint_fast64_t>(centLow*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t upper_ind = static_cast<uint_fast64_t>(centHigh*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t n_in_bin = upper_ind - lower_ind;

            //total number of events in this bin
            jet_file.write(reinterpret_cast<char*>(&n_in_bin), sizeof n_in_bin);

            auto it = colls_scatterings.crbegin();
            std::advance(it, lower_ind);

            for (uint_fast64_t ii = 0; ii<n_in_bin; it++, ii++)
            {
                for (auto e_co : it->second)
                {
                    auto e = e_co.dijet;

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                }
                io::append_single_coll_binary(jet_file, it->second, unirand, eng_shared);
            }
            jet_file.close();
            std::cout<<n_in_bin<<std::endl;

            dETdy_by_cent.add(new_ET_y);
            dEdy_by_cent.add(new_E_y);

            std::stringstream outname{""};

            outname<<"dEdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dEdy_by_cent, 
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
            outname.seekp(0);
            outname<<"dETdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dETdy_by_cent,
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
        }
    }

    if (save_endstate_jets)
    {
        std::vector<std::tuple<double,double> > centBins;
        std::ifstream cents_input(this->centrality_filename);
        if (!cents_input.is_open())
        {
            std::cout<<"Could not open "<<this->centrality_filename<<std::endl;
            centBins.push_back(std::tuple<double,double>{0.0, 1.0});
            return;
        }
        std::string line;
        for (; std::getline(cents_input, line);)
        {
            std::stringstream ss(line);
            std::string token;
            std::getline(ss, token, ',');
            double centLow = std::stod(token);
            std::getline(ss, token, ',');
            double centHigh = std::stod(token);
            centBins.push_back(std::tuple<double,double>{centLow, centHigh});
        }
        cents_input.close();

        std::string name_pfs{this->name+"_nomc.dat"};

        std::ofstream jet_file;

        for (auto [centLow, centHigh] : centBins)
        {
            std::stringstream jetsname{""};
            jetsname<<"jets_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            jet_file.open(jetsname.str(), std::ios::out | std::ios::binary);

            histo_1d dETdy_by_cent{y_bins};
            histo_1d dEdy_by_cent{y_bins};
            std::vector<std::tuple<double, double> > new_ET_y;
            std::vector<std::tuple<double, double> > new_E_y;

            uint_fast64_t N_evts_tot = colls_scatterings_nomc.size();
            double eps = 0.1/static_cast<double>(N_evts_tot);

            uint_fast64_t lower_ind = static_cast<uint_fast64_t>(centLow*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t upper_ind = static_cast<uint_fast64_t>(centHigh*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t n_in_bin = upper_ind - lower_ind;

            jet_file.write(reinterpret_cast<char*>(&n_in_bin), sizeof n_in_bin);

            auto it = colls_scatterings_nomc.crbegin();
            std::advance(it, lower_ind);

            for (uint_fast64_t ii = 0; ii<n_in_bin; it++, ii++)
            {
                for (auto e_co : it->second)
                {
                    auto e = e_co.dijet;

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                }
                io::append_single_coll_binary(jet_file, it->second, unirand, eng_shared);
            }
            jet_file.close();
            std::cout<<n_in_bin<<std::endl;

            dETdy_by_cent.add(new_ET_y);
            dEdy_by_cent.add(new_E_y);

            std::stringstream outname{""};

            outname<<"dEdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dEdy_by_cent, 
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
            outname.seekp(0);
            outname<<"dETdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dETdy_by_cent,
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
        }
    }

    if (save_endstate_jets)
    {
        std::vector<std::tuple<double,double> > centBins;
        std::ifstream cents_input(this->centrality_filename);
        if (!cents_input.is_open())
        {
            std::cout<<"Could not open "<<this->centrality_filename<<std::endl;
            centBins.push_back(std::tuple<double,double>{0.0, 1.0});
            return;
        }
        std::string line;
        for (; std::getline(cents_input, line);)
        {
            std::stringstream ss(line);
            std::string token;
            std::getline(ss, token, ',');
            double centLow = std::stod(token);
            std::getline(ss, token, ',');
            double centHigh = std::stod(token);
            centBins.push_back(std::tuple<double,double>{centLow, centHigh});
        }
        cents_input.close();

        std::string name_pfs{this->name+"_nosat.dat"};

        std::ofstream jet_file;

        for (auto [centLow, centHigh] : centBins)
        {
            std::stringstream jetsname{""};
            jetsname<<"jets_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            jet_file.open(jetsname.str(), std::ios::out | std::ios::binary);

            histo_1d dETdy_by_cent{y_bins};
            histo_1d dEdy_by_cent{y_bins};
            std::vector<std::tuple<double, double> > new_ET_y;
            std::vector<std::tuple<double, double> > new_E_y;

            uint_fast64_t N_evts_tot = colls_scatterings_nosat.size();
            double eps = 0.1/static_cast<double>(N_evts_tot);

            uint_fast64_t lower_ind = static_cast<uint_fast64_t>(centLow*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t upper_ind = static_cast<uint_fast64_t>(centHigh*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t n_in_bin = upper_ind - lower_ind;

            jet_file.write(reinterpret_cast<char*>(&n_in_bin), sizeof n_in_bin);

            auto it = colls_scatterings_nosat.crbegin();
            std::advance(it, lower_ind);

            for (uint_fast64_t ii = 0; ii<n_in_bin; it++, ii++)
            {
                for (auto e_co : it->second)
                {
                    auto e = e_co.dijet;

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                }
                io::append_single_coll_binary(jet_file, it->second, unirand, eng_shared);
            }
            jet_file.close();
            std::cout<<n_in_bin<<std::endl;

            dETdy_by_cent.add(new_ET_y);
            dEdy_by_cent.add(new_E_y);

            std::stringstream outname{""};

            outname<<"dEdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dEdy_by_cent, 
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
            outname.seekp(0);
            outname<<"dETdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dETdy_by_cent,
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
        }
    }

    if (save_endstate_jets)
    {
        std::vector<std::tuple<double,double> > centBins;
        std::ifstream cents_input(this->centrality_filename);
        if (!cents_input.is_open())
        {
            std::cout<<"Could not open "<<this->centrality_filename<<std::endl;
            centBins.push_back(std::tuple<double,double>{0.0, 1.0});
            return;
        }
        std::string line;
        for (; std::getline(cents_input, line);)
        {
            std::stringstream ss(line);
            std::string token;
            std::getline(ss, token, ',');
            double centLow = std::stod(token);
            std::getline(ss, token, ',');
            double centHigh = std::stod(token);
            centBins.push_back(std::tuple<double,double>{centLow, centHigh});
        }
        cents_input.close();

        std::string name_pfs{this->name+"_neither.dat"};

        std::ofstream jet_file;

        for (auto [centLow, centHigh] : centBins)
        {
            std::stringstream jetsname{""};
            jetsname<<"jets_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            jet_file.open(jetsname.str(), std::ios::out | std::ios::binary);

            histo_1d dETdy_by_cent{y_bins};
            histo_1d dEdy_by_cent{y_bins};
            std::vector<std::tuple<double, double> > new_ET_y;
            std::vector<std::tuple<double, double> > new_E_y;

            uint_fast64_t N_evts_tot = colls_scatterings_neither.size();
            double eps = 0.1/static_cast<double>(N_evts_tot);

            uint_fast64_t lower_ind = static_cast<uint_fast64_t>(centLow*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t upper_ind = static_cast<uint_fast64_t>(centHigh*static_cast<double>(N_evts_tot)+eps);
            uint_fast64_t n_in_bin = upper_ind - lower_ind;

            jet_file.write(reinterpret_cast<char*>(&n_in_bin), sizeof n_in_bin);

            auto it = colls_scatterings_neither.crbegin();
            std::advance(it, lower_ind);

            for (uint_fast64_t ii = 0; ii<n_in_bin; it++, ii++)
            {
                for (auto e_co : it->second)
                {
                    auto e = e_co.dijet;

                    new_ET_y.emplace_back(e.y1, e.kt);
                    new_ET_y.emplace_back(e.y2, e.kt);

                    new_E_y.emplace_back(e.y1, e.kt*cosh(e.y1));
                    new_E_y.emplace_back(e.y2, e.kt*cosh(e.y2));
                }
                io::append_single_coll_binary(jet_file, it->second, unirand, eng_shared);
            }
            jet_file.close();
            std::cout<<n_in_bin<<std::endl;

            dETdy_by_cent.add(new_ET_y);
            dEdy_by_cent.add(new_E_y);

            std::stringstream outname{""};

            outname<<"dEdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dEdy_by_cent, 
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
            outname.seekp(0);
            outname<<"dETdy_"<<static_cast<uint_fast16_t>(centLow*100)<<"_"<<static_cast<uint_fast16_t>(centHigh*100)<<"_"<<name_pfs;
            io::print_1d_histo
            (
                dETdy_by_cent,
                outname.str(), 
                1.0/ static_cast<double>(n_in_bin),
                false
            );
        }
    }
}

