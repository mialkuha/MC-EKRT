//Copyright (c) 2023 Mikko Kuha 

#include "mcaa.hpp"

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        mcaa sim("params_template");
        sim.run();  
    }
    else
    {
        const std::string filename{argv[1]};
        mcaa sim(filename);
        sim.run();
    }
    
    return 0;
}
