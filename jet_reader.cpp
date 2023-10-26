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
    int b1;
    int q1;
    int b2;
    int q2;
    dijet(const double &t01_, const double &t02_, const double &x_, const double &y_, 
            const double &pt_, const double &y1_, const double &y2_, const double &phi_, 
            const double &tata_, const int_fast8_t &f1, const int_fast8_t &f2)
    : t01(t01_), t02(t02_), x(x_), y(y_), pt(pt_), y1(y1_), y2(y2_), phi(phi_), tata(tata_)
    {
        if (f1 != 0)
        {
            b1 = static_cast<int>(f1 / abs(f1));
            if (f1 % 2 == 0)
            {
                q1 = 2*b1;
            }
            else
            {
                q1 = -1*b1;
            }
        }
        else
        {
            b1 = 0;
            q1 = 0;
        }
        if (f2 != 0)
        {
            b2 = static_cast<int>(f2 / abs(f2));
            if (f2 % 2 == 0)
            {
                q2 = 2*b2;
            }
            else
            {
                q2 = -1*b2;
            }
        }
        else
        {
            b2 = 0;
            q2 = 0;
        }
    }
};

 
int main()
{
    std::string filename = "jets_0_100_example_name_new_nool.dat";
    std::ifstream jet_file(filename, std::ios::in | std::ios::binary);

    if (!jet_file.is_open())
    {
        std::cout << "Error opening file " << filename << std::endl;
        exit(1);
    }

    uint_fast64_t n_events, n_jets;
    double t01, t02, x, y, pt, y1, y2, phi, tata;
    int_fast8_t f1, f2;
    int quarks, charges;
    double total_B = 0;
    double total_Q = 0;
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

    uint_fast16_t n_bins = 50;

    uint_fast64_t n_dijet_total = 0;
    for (uint_fast64_t i=0; i<n_events; i++)
    {
        jet_file.read(reinterpret_cast<char*>(&n_jets), sizeof n_jets);

        std::vector<dijet> jets;
        jets.reserve(n_jets);
        quarks = 0;
        charges = 0;

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
            jet_file.read(reinterpret_cast<char*>(&f1) , sizeof f1); 
            jet_file.read(reinterpret_cast<char*>(&f2) , sizeof f2); 

            jets.emplace_back(t01, t02, x, y, pt, y1, y2, phi, tata, f1, f2);
        }
        events_jets.emplace_back(std::move(jets));
    }
    jet_file.close();

    std::cout << "Done! Read " << events_jets.size() << " events totaling "
              << n_dijet_total << " dijets, average baryon number was "
              << total_B << " and average em-charge was "
              << total_Q <<std::endl;

    double delta = (0 + 11.0)/n_bins;
    for (int i = 0; i < n_bins+1; i++)
    {
        std::cout << i*delta << ',';
    }
    std::cout << std::endl;

    std::vector<double> tata_in_bin(n_bins, 0.0);
    std::vector<int> jets_in_bin(n_bins, 0);

    double tata_av = 0.0;

    for (auto event : events_jets)
    {
        for (auto dijet : event)
        {
            bool found =false;
            double x = dijet.x;
            double y = dijet.y;
            double r = std::sqrt(x*x + y*y);
            for (int i = 0; i < n_bins-1; i++)
            {
                if ( (i+1)*delta > r)
                {
                    found = true;
                    tata_in_bin[i] += dijet.tata;
                    jets_in_bin[i]++;
                    break;
                }
            }
            
            if (!found)
            {
                tata_in_bin[n_bins-1] += dijet.tata;
                jets_in_bin[n_bins-1]++;
            }
        }
    }
        
    for (int i = 0; i < n_bins; i++)
    {
        std::cout << tata_in_bin[i] / static_cast<double>(jets_in_bin[i]* /*events_jets.size()**/delta/**((i+0.5)*delta)*/) << ',';
    }
    std::cout << std::endl;
}