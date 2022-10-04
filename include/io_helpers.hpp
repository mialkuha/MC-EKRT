//Copyright (c) 2022 Mikko Kuha

#ifndef IO_HELPERS_HPP
#define IO_HELPERS_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "histo.hpp"
#include "nucleon.hpp"
#include "typedefs.hpp"

class io
{
public:
    class Coll 
    {

        private:
        double b, x, nTwo, sumET;
        ulong nColl, nPart, nEnd;
        bool alice = true;

        public:
        Coll() { b = 0.; x = 0.; nTwo = 0.; nColl = 0; nPart = 0; nEnd = 0; sumET = 0.;}
        Coll(ulong nColl_, ulong nPart_, ulong nEnd_, double b_, double sumET_){
            b = b_; nColl = nColl_; nPart = nPart_; nEnd = nEnd_; sumET = sumET_;
            // ATLAS definitions. 
            if (!alice) {
            x = 0.09;
            nTwo = 0.5 * (1. - x) * static_cast<double>(nPart) + x * static_cast<double>(nColl);
            // ALICE definitions.
        } else {
            x = 0.801;
            nTwo = x * static_cast<double>(nPart) + (1. - x) * static_cast<double>(nColl);
        }
        }
        Coll(ulong nCollIn, ulong nPartIn, ulong nEndIn, double bIn, std::shared_ptr<std::mt19937> random_generator){
            b = bIn; nColl = nCollIn; nPart = nPartIn; nEnd = nEndIn; sumET = 0.;
            // ATLAS definitions. 
            if (!alice) {
            x = 0.09; 
            nTwo = 0.5 * (1. - x) * static_cast<double>(nPart) + x * static_cast<double>(nColl); sampleET(random_generator);
            // ALICE definitions.
        } else {
            x = 0.801;
            nTwo = x * static_cast<double>(nPart) + (1. - x) * static_cast<double>(nColl); sampleET(random_generator);
        }
        }
        double getB() const {return b;}
        ulong getNcoll() const {return nColl;}
        ulong getNpart() const {return nPart;}
        ulong getNend() const {return nEnd;}
        double getNtwo() const {return nTwo;}
        double getET() const {return sumET;}
        void sampleET(std::shared_ptr<std::mt19937> random_generator) {
            double mu = 46.4;
            double k = 1.5;
            double p = (mu/k) / (mu/k + 1.);
            uint nAncestors = static_cast<uint>(floor(nTwo));
            uint mult = 0;
            //std::cout << (*random_generator)() << ' ' << std::flush;
            std::negative_binomial_distribution<uint> distr(static_cast<uint>(k),1.-p);
            for (uint i = 0; i < nAncestors; ++i) {
            mult += distr(*random_generator);
            //std::cout << mult << ' ' << std::flush;
            }
            sumET = mult;
            //std::cout << std::endl << nAncestors << ' ' << sumET/nAncestors << std::endl << std::endl;
            // cout << endl << nAncestors << ' ' << sumET/nAncestors << endl << endl;
        }
    };
    static bool compNpart(Coll c1, Coll c2){
        return c1.getNpart() < c2.getNpart();
        }
    static bool compNcoll(Coll c1, Coll c2){
        return c1.getNcoll() < c2.getNcoll();
        }
    static bool compNend(Coll c1, Coll c2){
        return c1.getNend() < c2.getNend();
        }
    static bool compNtwo(Coll c1, Coll c2){
        return c1.getNtwo() < c2.getNtwo();
        }
    static bool compET(Coll c1, Coll c2){
        return c1.getET() < c2.getET();
        }
    template <typename T>
    static double calc_ave(std::vector<Coll> &collisions, T (Coll::*func)()const ){
        std::vector<Coll>::iterator it;
        T tot = 0;
        for (it = collisions.begin(); it != collisions.end(); ++it){
            tot += ((*it).*func)();
        }
        return static_cast<double>(tot)/static_cast<double>(collisions.size());  
    }


    static auto mc_glauber_style_report
    (
        std::vector<Coll> &collisions, 
        const xsectval &sigma_inel, 
        const uint &N_events, 
        const uint &nBins, 
        const double *const binsLow, 
        const double *const binsHigh,
        std::ostream &out_stream
    ) noexcept -> void
    {
        // Make sure that no rounding downwards.
        double eps = 0.1/N_events;
        
        // int comp = 1;
        for (uint8_t comp = 1; comp <= 5; ++comp) {
            
            // Print out the header.
            out_stream << std::endl;
            out_stream << "# sigma_inel = " << sigma_inel << " mb" << std::endl;
            
            switch (comp)
            {
                case 1:
                    std::sort(collisions.begin(), collisions.end(), compNpart);
                    out_stream << "Using Npart for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     nPmax     nPmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                case 2:
                    std::sort(collisions.begin(), collisions.end(), compNcoll);
                    out_stream << "Using Ncoll for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     nCmax     nCmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                case 3:
                    std::sort(collisions.begin(), collisions.end(), compNtwo);
                    out_stream << "Using two-component (ancestors) model for centrality determination"
                            << std::endl << std::endl
                            << "#  cent     nAmax     nAmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                case 4:
                    std::sort(collisions.begin(), collisions.end(), compET);
                    out_stream << "Using sumET model for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     ETmax     ETmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
                    break;
                default: //case 5
                    std::sort(collisions.begin(), collisions.end(), compET);
                    out_stream << "Using the number of produced jets for centrality determination" 
                            << std::endl << std::endl
                            << "#  cent     nJmax     nJmin       <b>   <nPart>   <nColl>    <nEnd>    <T_AA>     <E_T>"
                            << std::endl;
            }
            std::reverse(collisions.begin(), collisions.end());
            
            std::vector<Coll> centrality;
            // Number of collisions within centrality class i.
            for (uint i = 0; i < nBins; i++)
            {
                uint lower = static_cast<uint>(binsLow[i]*N_events+eps);
                uint upper = static_cast<uint>(binsHigh[i]*N_events+eps);
                centrality.clear();
                while (lower<upper)
                {
                    centrality.push_back(collisions.at(lower++));
                }
                
                // Print the centrality selection.
                out_stream << std::setw(2) << binsLow[i]*100 << " " << std::setw(4) << binsHigh[i]*100;
                
                switch (comp)
                {
                    case 1:
                        out_stream << std::setw(10) << centrality[0].getNpart()
                                   << std::setw(10) << centrality[centrality.size() - 1].getNpart();
                        break;
                    case 2:
                        out_stream << std::setw(10) << centrality[0].getNcoll() 
                                   << std::setw(10) << centrality[centrality.size() - 1].getNcoll();
                        break;
                    case 3:
                        out_stream << std::setw(10) <<centrality[0].getNtwo()
                                   << std::setw(10) << centrality[centrality.size() - 1].getNtwo();
                        break;
                    case 4:
                        out_stream << std::setw(10) << centrality[0].getET() 
                                   << std::setw(10) << centrality[centrality.size() - 1].getET();
                        break;
                    default : //case 5
                        out_stream << std::setw(10) << centrality[0].getNend() 
                                   << std::setw(10) << centrality[centrality.size() - 1].getNend();
                        break;
                }
                
                out_stream << std::setw(10) << io::calc_ave<double>(centrality, &Coll::getB) 
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNpart) 
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNcoll)
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNend) 
                           << std::setw(10) << io::calc_ave<ulong>(centrality, &Coll::getNcoll)/sigma_inel
                           << std::setw(10) << calc_ave<double>(centrality, &Coll::getET)
                           << std::endl;      
            }
        }
    }

    static auto print_3d_histo
    (
        const histo_3d &histo_, 
        const std::string &filename,
        const double normalization = 1.0,
        const bool divide_tot = true
    ) noexcept -> void
    {
        std::ofstream file;

        file.open(filename);

        auto [ histo, uf, of, xs, tot ] = histo_.project_2d(1);
        auto [ xs0, xs1 ] = xs;

        file << "///Total count: "<<tot<<std::endl;
        file << "///underflow: "<<uf<<std::endl;
        file << "///overflow: "<<of<<std::endl;

        file << "///y bin walls:  ";
        for (auto y : xs1)
        {
            file<<y<<' ';
        }
        file << std::endl;

        auto norm = (divide_tot)? normalization : normalization * static_cast<double>(tot);

        ulong n_x_bins = xs0.size()-1;
        ulong n_y_bins = xs1.size()-1;

        for (ulong i = 0; i < n_x_bins; i++)
        {
            file << xs0[i] << ' ' << xs0[i+1] << ' ';

            for (ulong j = 0; j < n_y_bins; j++)
            {   
                file << histo[i][j] * norm << ' ';
            }
            file << std::endl;
        }

        auto [ histo2, uf2, of2, xs2, tot2 ] = histo_.project_2d(2);
        auto [ xs02, xs12 ] = xs2;
        file << "///y bin walls:  ";
        for (auto y : xs12)
        {
            file<<y<<' ';
        }
        file << std::endl;

        norm = (divide_tot)? normalization : normalization * static_cast<double>(tot2);

        n_x_bins = xs02.size()-1;
        n_y_bins = xs12.size()-1;

        for (ulong i = 0; i < n_x_bins; i++)
        {
            file << xs02[i] << ' ' << xs02[i+1] << ' ';

            for (ulong j = 0; j < n_y_bins; j++)
            {   
                file << histo2[i][j] * norm << ' ';
            }
            file << std::endl;
        }
        std::cout<<"printed to "<<filename<<std::endl;
        file.close();
    }

    static auto print_2d_histo
    (
        const histo_2d &histo_, 
        const std::string &filename,
        const double normalization = 1.0,
        const bool divide_tot = true
    ) noexcept -> void
    {
        std::ofstream file;

        file.open(filename);

        auto [ histo, uf, of, xs, tot ] = histo_.get_histo();
        auto [ uf0, uf1 ] = uf;
        auto [ of0, of1 ] = of;
        auto [ xs0, xs1 ] = xs;

        file << "///Total count: "<<tot<<std::endl;
        file << "///underflow: "<<uf0<<' '<<uf1<<std::endl;
        file << "///overflow: "<<of0<<' '<<of1<<std::endl;
        file << "///normalization: "<<normalization<<std::endl;

        file << "///y bin walls:  ";
        for (auto y : xs1)
        {
            file<<y<<' ';
        }
        file << std::endl;

        auto norm = (divide_tot)? normalization : normalization * static_cast<double>(tot);

        auto n_x_bins = histo.size();
        auto n_y_bins = histo[0].size();

        for (uint8_t i = 0; i < n_x_bins; i++)
        {
            file << xs0[i] << ' ' << xs0[i+1] << ' ';

            for (uint8_t j = 0; j < n_y_bins; j++)
            {   
                file << histo[i][j] * norm << ' ';
            }
            file << std::endl;
        }
        std::cout<<"printed to "<<filename<<std::endl;
        file.close();
    }

    static auto save_event
    (
        std::ofstream &event_file, 
        const std::vector<nucleon> &pro, 
        const std::vector<nucleon> &tar, 
        const spatial &impact_parameter
    ) noexcept -> void
    {
        event_file << "{"<<std::endl;
        event_file << "    "<<impact_parameter<<','<<std::endl;
        event_file << "    {"<<std::endl;
        for (const auto & n : pro)
        {
            event_file << "        {"<<std::endl;
            event_file << "            "<<n.co<<','<<std::endl;
            event_file << "            "<<n.mom<<','<<std::endl;
            event_file << "            "<<n.wounded<<','<<std::endl;
            event_file << "            "<<n.is_neutron<<std::endl;
            if (&n != &pro.back())
            {
                event_file << "        },"<<std::endl;
            }
            else
            {
                event_file << "        }"<<std::endl;
            }
        }
        event_file << "    },"<<std::endl;
        event_file << "    {"<<std::endl;
        for (auto n : tar)
        {
            event_file << "        {"<<std::endl;
            event_file << "            "<<n.co<<std::endl;
            event_file << "            "<<n.mom<<std::endl;
            event_file << "            "<<n.wounded<<std::endl;
            event_file << "            "<<n.is_neutron<<std::endl;
            if (&n != &tar.back())
            {
                event_file << "        },"<<std::endl;
            }
            else
            {
                event_file << "        }"<<std::endl;
            }
        }
        event_file << "    }"<<std::endl;
        event_file << "},"<<std::endl;
    }

    static auto save_event
    (
        std::ofstream &event_file, 
        const std::vector<nucleon> &pro, 
        const std::vector<nucleon> &tar, 
        const spatial &impact_parameter,
        const std::vector<dijet_specs> &filtered_scatterings
    ) noexcept -> void
    {
        event_file << "{"<<std::endl;
        event_file << "    "<<impact_parameter<<','<<std::endl;
        event_file << "    {"<<std::endl;
        for (const auto & n : pro)
        {
            event_file << "        {"<<std::endl;
            event_file << "            "<<n.co<<','<<std::endl;
            event_file << "            "<<n.mom<<','<<std::endl;
            event_file << "            "<<n.wounded<<','<<std::endl;
            event_file << "            "<<n.is_neutron<<std::endl;
            if (&n != &pro.back())
            {
                event_file << "        },"<<std::endl;
            }
            else
            {
                event_file << "        }"<<std::endl;
            }
        }
        event_file << "    },"<<std::endl;
        event_file << "    {"<<std::endl;
        for (const auto & n : tar)
        {
            event_file << "        {"<<std::endl;
            event_file << "            "<<n.co<<','<<std::endl;
            event_file << "            "<<n.mom<<','<<std::endl;
            event_file << "            "<<n.wounded<<','<<std::endl;
            event_file << "            "<<n.is_neutron<<std::endl;
            if (&n != &tar.back())
            {
                event_file << "        },"<<std::endl;
            }
            else
            {
                event_file << "        }"<<std::endl;
            }
        }
        event_file << "    },"<<std::endl;
        event_file << "    {"<<std::endl;
        for (const auto & sc : filtered_scatterings)
        {
            event_file << "        {"<<std::endl;
            event_file << "            "<<sc.kt<<','<<std::endl;
            event_file << "            "<<sc.y1<<','<<std::endl;
            event_file << "            "<<sc.y2<<','<<std::endl;
            event_file << "            "<<sc.init1<<','<<std::endl;
            event_file << "            "<<sc.init2<<','<<std::endl;
            event_file << "            "<<sc.final1<<','<<std::endl;
            event_file << "            "<<sc.final2<<std::endl;
            if (&sc != &filtered_scatterings.back())
            {
                event_file << "        },"<<std::endl;
            }
            else
            {
                event_file << "        }"<<std::endl;
            }
        }
        event_file << "    }"<<std::endl;
        event_file << "},"<<std::endl;
    }

    static auto print_1d_histo
    (
        const histo_1d &histo_, 
        const std::string &filename,
        const double normalization = 1.0,
        const bool divide_tot = true
    ) noexcept -> void
    {
        std::ofstream file;

        file.open(filename);

        auto [ histo, uf, of, xs, tot ] = histo_.get_histo();

        file << "///Total count: "<<tot<<std::endl;
        file << "///underflow: "<<uf<<std::endl;
        file << "///overflow: "<<of<<std::endl;

        file << "///x bin walls:  ";
        for (auto x : xs)
        {
            file<<x<<' ';
        }
        file << std::endl;

        auto norm = (divide_tot)? normalization : normalization * static_cast<double>(tot);

        for (auto e : histo)
        {
            file<< e * norm <<' ';
        }
        file << std::endl;

        std::cout<<"printed to "<<filename<<std::endl;
        file.close();
    }

    static auto read_conf
    (
        const std::string filename
    ) noexcept
    {
        std::ifstream input(filename);
        
        std::string name;
        bool is_pp              = false;
        bool is_pa              = false;
        bool is_aa              = true;
        bool use_npdfs          = false;
        bool use_snpdfs         = false;
        bool is_mom_cons        = false;
        bool is_mom_cons_new    = false;
        bool are_ns_depleted    = false;
        bool is_saturation      = false;
        bool is_sat_overlap     = false;
        uint32_t n_events       = 10000;
        spatial b_max           = 20;//fm
        momentum pt0            = 2.728321;//GeV

        uint16_t count = 0;

        if (input.is_open())
        {
            std::string line;

            for (; std::getline(input, line); )
            {
                std::istringstream line_stream(line);
                std::string param_name;
                line_stream >> param_name;

                if (param_name == "name")
                {
                    line_stream >> name;
                }
                else if (param_name == "is_pp")
                {
                    line_stream >> std::boolalpha >> is_pp;
                }
                else if (param_name == "is_pa")
                {
                    line_stream >> std::boolalpha >> is_pa;
                }
                else if (param_name == "is_aa")
                {
                    line_stream >> std::boolalpha >> is_aa;
                }
                else if (param_name == "use_npdfs")
                {
                    line_stream >> std::boolalpha >> use_npdfs;
                }
                else if (param_name == "use_snpdfs")
                {
                    line_stream >> std::boolalpha >> use_snpdfs;
                }
                else if (param_name == "is_mom_cons")
                {
                    line_stream >> std::boolalpha >> is_mom_cons;
                }
                else if (param_name == "is_mom_cons_new")
                {
                    line_stream >> std::boolalpha >> is_mom_cons_new;
                }
                else if (param_name == "are_ns_depleted")
                {
                    line_stream >> std::boolalpha >> are_ns_depleted;
                }
                else if (param_name == "is_saturation")
                {
                    line_stream >> std::boolalpha >> is_saturation;
                }
                else if (param_name == "is_sat_overlap")
                {
                    line_stream >> std::boolalpha >> is_sat_overlap;
                }
                else if (param_name == "n_events")
                {
                    line_stream >> n_events;
                }
                else if (param_name == "b_max")
                {
                    line_stream >> b_max;
                }
                else if (param_name == "pt0")
                {
                    line_stream >> pt0;
                }
                count++;
            }
        }
        else
        {
            std::cout<<"Error opening the configuration file!!"<<std::endl;
            std::terminate();
        }

        std::cout<<"Read "<< count <<" parameters from a file "<< filename <<std::endl;

        return std::make_tuple
                (
                    name,
                    is_pp,
                    is_pa,
                    is_aa,
                    use_npdfs,
                    use_snpdfs,
                    is_mom_cons,
                    is_mom_cons_new,
                    are_ns_depleted,
                    is_saturation,
                    is_sat_overlap,
                    n_events,
                    b_max,
                    pt0
                );
    }

    static auto print_histos
    (
        const std::string &name_postfix,
        const std::vector<histo_2d> &jets,
        const std::vector<histo_2d> &dijets,
        const std::vector<histo_1d> &dETdy,
        const std::vector<histo_1d> &dEdy,
        const std::vector<histo_1d> &dNdy,
        const std::vector<histo_1d> &dNdET,
        const std::vector<histo_1d> &dETdeta,
        const std::vector<histo_1d> &dEdeta,
        const std::vector<histo_1d> &dETdb,
        const std::vector<histo_1d> &dEdb,
        const xsectval &dijet_norm,
        const uint32_t &AA_events
    ) noexcept -> void
    {
        std::array<std::string, 6> true_postfixes{name_postfix+".dat",
                                                  name_postfix+"_MC.dat",
                                                  name_postfix+"_MC_ND.dat",
                                                  name_postfix+"_SAT.dat",
                                                  name_postfix+"_SAT_OL.dat",
                                                  name_postfix+"_SAT_MC.dat"};
        
        for (size_t i=0; i<6; i++)
        { //sigma1jet
            print_2d_histo
            (
                jets[i], 
                "sigma1jet_sim_"+true_postfixes[i], 
                2.0 * dijet_norm
            );
        }

        for (size_t i=0; i<6; i++)
        { //dNdpTdy
            print_2d_histo
            (
                jets[i],
                "dNdpTdy_sim_"+true_postfixes[i], 
                1.0,
                false
            );
        }

        for (size_t i=0; i<6; i++)
        { //sigmadijet
            print_2d_histo
            (
                dijets[i],
                "sigmadijet_sim_"+true_postfixes[i], 
                dijet_norm
            );
        }

        for (size_t i=0; i<6; i++)
        { //dETdy
            print_1d_histo
            (
                dETdy[i],
                "dETdy_sim_"+true_postfixes[i], 
                1.0/ AA_events,
                false
            );
        }

        for (size_t i=0; i<6; i++)
        { //dEdy
            print_1d_histo
            (
                dEdy[i], 
                "dEdy_sim_"+name_postfix+true_postfixes[i], 
                1.0/ AA_events,
                false
            );
        }

        for (size_t i=0; i<6; i++)
        { //dNdy
            print_1d_histo
            (
                dNdy[i], 
                "dNdy_sim_"+true_postfixes[i], 
                1.0/ AA_events,
                false
            );
        }

        for (size_t i=0; i<6; i++)
        { //dNdET
            print_1d_histo
            (
                dNdET[i], 
                "dNdET_sim_"+true_postfixes[i], 
                1.0,
                false
            );
        }
        
        for (size_t i=0; i<6; i++)
        { //dETdeta
            print_1d_histo
            (
                dETdeta[i], 
                "dETdeta_sim_"+true_postfixes[i], 
                1.0 / AA_events,
                false
            );
        }

        for (size_t i=0; i<6; i++)
        { //dEdeta
            print_1d_histo
            (
                dEdeta[i], 
                "dEdeta_sim_"+true_postfixes[i], 
                1.0 / AA_events,
                false
            );
        }

        for (size_t i=0; i<6; i++)
        { //dETdb
            print_1d_histo
            (
                dETdb[i], 
                "dETdb_sim_"+true_postfixes[i], 
                1.0,
                true
            );
        }

        for (size_t i=0; i<6; i++)
        { //dEdb
            print_1d_histo
            (
                dEdb[i],
                "dEdb_sim_"+true_postfixes[i], 
                1.0,
                true
            );
        }
    }

    static auto print_histos
    (
        const std::string &name_postfix,
        const histo_2d &jets,
        const histo_2d &dijets,
        const histo_1d &dETdy,
        const histo_1d &dEdy,
        const histo_1d &dNdy,
        const histo_1d &dNdET,
        const histo_1d &dETdeta,
        const histo_1d &dEdeta,
        const histo_1d &dETdb,
        const histo_1d &dEdb,
        const xsectval &dijet_norm,
        const uint32_t &AA_events
    ) noexcept -> void
    {
        std::string true_postfixes{name_postfix+".dat"};
        
        print_2d_histo
        (
            jets, 
            "sigma1jet_sim_"+true_postfixes, 
            2.0 * dijet_norm
        );
        
        print_2d_histo
        (
            jets,
            "dNdpTdy_sim_"+true_postfixes, 
            1.0,
            false
        );
        
        print_2d_histo
        (
            dijets,
            "sigmadijet_sim_"+true_postfixes, 
            dijet_norm
        );
        
        print_1d_histo
        (
            dETdy,
            "dETdy_sim_"+true_postfixes, 
            1.0/ AA_events,
            false
        );
        
        print_1d_histo
        (
            dEdy, 
            "dEdy_sim_"+name_postfix+true_postfixes, 
            1.0/ AA_events,
            false
        );
        
        print_1d_histo
        (
            dNdy, 
            "dNdy_sim_"+true_postfixes, 
            1.0/ AA_events,
            false
        );
        
        print_1d_histo
        (
            dNdET, 
            "dNdET_sim_"+true_postfixes, 
            1.0,
            false
        );
        
        print_1d_histo
        (
            dETdeta, 
            "dETdeta_sim_"+true_postfixes, 
            1.0 / AA_events,
            false
        );
        
        print_1d_histo
        (
            dEdeta, 
            "dEdeta_sim_"+true_postfixes, 
            1.0 / AA_events,
            false
        );
        
        print_1d_histo
        (
            dETdb, 
            "dETdb_sim_"+true_postfixes, 
            1.0,
            true
        );
        
        print_1d_histo
        (
            dEdb,
            "dEdb_sim_"+true_postfixes, 
            1.0,
            true
        );
        
    }
};

#endif // IO_HELPERS_HPP