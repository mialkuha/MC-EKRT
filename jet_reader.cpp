#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

struct dijet
{
    double t01;
    double t02;
    double x;
    double y;
    double pt;
    double y1;
    double y2;
    double phi;
    double tata;
    dijet(const double &t01_, const double &t02_, const double &x_, const double &y_, 
            const double &pt_, const double &y1_, const double &y2_, const double &phi_, 
            const double &tata_)
    : t01(t01_), t02(t02_), x(x_), y(y_), pt(pt_), y1(y1_), y2(y2_), phi(phi_) , tata(tata_) { }
};

 
int main()
{
    std::string filename = "jets_0_5_LHC_onoff.dat";
    std::ifstream jet_file(filename, std::ios::in | std::ios::binary);

    if (!jet_file.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    uint_fast64_t n_events, n_jets;
    double t01, t02, x, y, pt, y1, y2, phi, tata;
    std::vector<std::vector<dijet> > events_jets;
    
    if (sizeof n_jets != 8 || sizeof t01 != 8 )
    {
        std::cout << "Change the types! n_jets should be 64 bit unsigned integer"
                  << " and the coordinates should all be 64 bit float numbers" << std::endl;
        exit(1);
    }

    jet_file.read(reinterpret_cast<char*>(&n_events), sizeof n_events);
    std::cout << "Trying to read " << n_events << " events from the file " 
              << filename << " ..." << std::endl;
    events_jets.reserve(n_events);

    uint_fast16_t n_bins = 20;

    uint_fast64_t n_dijet_total = 0;
    for (uint_fast64_t i=0; i<n_events; i++)
    {
        jet_file.read(reinterpret_cast<char*>(&n_jets), sizeof n_jets);

        std::vector<dijet> jets;
        jets.reserve(n_jets);

        for (uint_fast64_t ii=0; ii<n_jets; ii++, n_dijet_total++)
        {
            jet_file.read(reinterpret_cast<char*>(&t01)  , sizeof t01);  
            jet_file.read(reinterpret_cast<char*>(&t02)  , sizeof t02);  
            jet_file.read(reinterpret_cast<char*>(&x)    , sizeof x);    
            jet_file.read(reinterpret_cast<char*>(&y)    , sizeof y);    
            jet_file.read(reinterpret_cast<char*>(&pt)   , sizeof pt);   
            jet_file.read(reinterpret_cast<char*>(&y1)   , sizeof y1);   
            jet_file.read(reinterpret_cast<char*>(&y2)   , sizeof y2);   
            jet_file.read(reinterpret_cast<char*>(&phi)  , sizeof phi);  
            jet_file.read(reinterpret_cast<char*>(&tata) , sizeof tata); 

            jets.emplace_back(t01, t02, x, y, pt, y1, y2, phi, tata);
        }
        events_jets.emplace_back(std::move(jets));
    }
    jet_file.close();

    std::cout << "Done! Read " << events_jets.size() << " events totaling "
              << n_dijet_total << " dijets" <<std::endl;

    double delta = (20.0-1.0)/n_bins;
    std::ofstream dNdpt;
    dNdpt.open("dNdpt.csv");
    for (int i = 0; i < n_bins+1; i++)
    {
        dNdpt << 1.0+i*delta << ',';
    }
    dNdpt << std::endl;

    std::vector<std::vector<double> > pts_in_bin
    (
        n_bins+1, 
        std::vector<double>()
    );

    for (auto event : events_jets)
    {
        for (auto dijet : event)
        {
            bool found =false;
            for (int i = 0; i < n_bins; i++)
            {
                if ( 1.0+(i+1)*delta > dijet.pt)
                {
                    found = true;
                    pts_in_bin[i].push_back(dijet.pt);
                    break;
                }
            }
            if (!found)
            {
                pts_in_bin[n_bins].push_back(dijet.pt);
            }
        }
    }
        
    for (int i = 0; i < n_bins+1; i++)
    {
        dNdpt << pts_in_bin[i].size() << ',';
    }
    dNdpt << std::endl;
    dNdpt.close();
}

