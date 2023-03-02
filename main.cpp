//Copyright (c) 2022 Mikko Kuha 
//TODO get rid of all the hardcoded "208"s

#include "mcaa.hpp"

int main()
{
    const std::array<double, 5> K_s{4.0,8.0};
    const std::array<double, 5> M_s{1.0,2.0,4.0,8.0,16.0};
    for (auto K : K_s)
    {
        for (auto M : M_s)
        {
            std::cout<<std::endl<<std::endl<<K<<' '<<M<<':'<<std::endl<<std::endl<<std::endl;
            mcaa sim("params_template");
            sim.K_factor = K;
            sim.M_factor = M;
            sim.diff_params.K_factor = K;
            sim.jet_params.d_params.K_factor = K;
            sim.run();
        }
    }
    
    return 0;
}