#include "problem7_functions.hpp"
#include "problem9_functions.hpp"
#include "problem10_functions.hpp"
#include <iostream>
#include <chrono>

double timer_general(int n_step){
    auto tstart = std::chrono::high_resolution_clock::now();
    thomas_algo(n_step, -1, 2, -1);
    auto tstop = std::chrono::high_resolution_clock::now();
    double thomas_general_duration_seconds = std::chrono::duration<double>(tstop - tstart).count();
    return thomas_general_duration_seconds;
}

double timer_special(int n_step){
    auto tstart = std::chrono::high_resolution_clock::now();
    special_thomas_algo(n_step);
    auto tstop = std::chrono::high_resolution_clock::now();
    double thomas_special_duration_seconds = std::chrono::duration<double>(tstop - tstart).count();
    return thomas_special_duration_seconds;
}