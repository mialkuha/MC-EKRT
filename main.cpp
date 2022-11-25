//Copyright (c) 2022 Mikko Kuha 
//TODO implement NN-tree 

#include "mcaa.hpp"

int main()
{
    std::ofstream paramfile;
    {
    std::string name{"test_sigma_coupled6"};
    paramfile.open(name);
    paramfile<<"name "<<name<<std::endl;
    paramfile<<"sigma_inel_from_sigma_jet true"<<std::endl;
    paramfile<<"is_saturation false"<<std::endl;
    paramfile<<"is_mom_cons false"<<std::endl;
    paramfile<<"K_factor 2.0"<<std::endl;
    paramfile<<"K_factor 2.0"<<std::endl;
    paramfile<<"n_events 5000"<<std::endl;
    paramfile.close();

    mcaa sim(name);
    sim.run();
    }
    {
    std::string name{"test_sigma_uncoupled6"};
    paramfile.open(name);
    paramfile<<"name "<<name<<std::endl;
    paramfile<<"sigma_inel_from_sigma_jet false"<<std::endl;
    paramfile<<"is_saturation false"<<std::endl;
    paramfile<<"is_mom_cons false"<<std::endl;
    paramfile<<"K_factor 2.0"<<std::endl;
    paramfile<<"M_factor 5.0"<<std::endl;
    paramfile<<"n_events 5000"<<std::endl;
    paramfile.close();

    mcaa sim(name);
    sim.run();
    }
    return 0;
}