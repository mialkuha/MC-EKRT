//Copyright (c) 2021 Mikko Kuha

#ifndef PQCD_HPP
#define PQCD_HPP

#include "nn_coll.hpp"
#include "typedefs.hpp"

class nn_coll;//Declaration to get away from the circular references
class pqcd
{
public:
    static void generate_bin_NN_coll(nn_coll * coll) noexcept;
protected:
private:
};

#endif // PQCD_HPP
