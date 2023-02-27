//Copyright (c) 2022 Mikko Kuha 
//TODO get rid of all the hardcoded "208"s

#include "mcaa.hpp"

int main()
{
    mcaa sim("params_template");
    sim.run();
    return 0;
}