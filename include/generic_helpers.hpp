//Copyright (c) 2022 Mikko Kuha

#ifndef GENERIC_HELPERS_HPP
#define GENERIC_HELPERS_HPP

#include <limits>
#include <vector>

#include "typedefs.hpp"

class helpers
{
public:
    // return an evenly spaced 1d grid of doubles.
    static auto linspace
    (
        const double &first, 
        const double &last, 
        const uint16_t &len
    ) noexcept -> std::vector<double>
    {
        std::vector<double> result(len);
        double step = (last-first) / (len - 1);
        for (uint16_t i=0; i<len; i++) { result[i] = first + i*step; }
        return result;
    }

    // return an  1d grid of doubles distributed evenly in logarithm
    static auto loglinspace
    (
        const double &first, 
        const double &last, 
        const uint16_t &len
    ) noexcept -> std::vector<double>
    {
        std::vector<double> result(len);
        double step = pow(last/first, 1.0/(len-1.0));
        result[0] = first;
        for (uint16_t i=1; i<len; i++) { result[i] = result[i-1]*step; }
        return result;
    }

    //Finds a zero of the function f
    template<typename F, typename Ret_type, typename Arg_type>
    static auto secant_method
    (
        Arg_type *const x, 
        F f, 
        const Ret_type error_tolerance, 
        Ret_type *const last_fx,
        const Arg_type lambda=1.0 
    ) -> void
    {
        Arg_type x_n = *x, x_np1 = 1.2*x_n, x_np2 = 0.0;
        Ret_type fx_n = f(x_n);

        if (abs(fx_n) < error_tolerance)
        {
            *last_fx = fx_n;
            return;
        }

        Ret_type fx_np1 = f(x_np1);

        while ( abs(fx_np1) > error_tolerance )
        {
            x_np2 = x_np1 - lambda*fx_np1*(x_np1 - x_n)/(fx_np1 - fx_n);

            if (x_np2 <= 0.0) x_np2=1.0/std::numeric_limits<Arg_type>::max();

            fx_n = fx_np1;
            x_n = x_np1;
            x_np1 = x_np2;

            fx_np1 = f(x_np1);

            if (x_n == x_np1)
            {
                std::cout<<"Doesn't converge!!!"<<std::endl;
                std::cout<<x_n<<' '<<x_np1<<' '<<fx_np1<<std::endl;
                return;
            }
        }

        *last_fx = fx_np1;
        *x = x_np1;
        return;
    }


};

#endif // GENERIC_HELPERS_HPP