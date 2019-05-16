/**
*      _________   _____________________  ____  ______
*     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
*    / /_  / /| | \__ \ / / / /   / / / / / / / __/
*   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
*  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
*
*  http://www.acl.inf.ethz.ch/teaching/fastcode
*  How to Write Fast Numerical Code 263-2300 - ETH Zurich
*  Copyright (C) 2019 
*                   Tyler Smith        (smitht@inf.ethz.ch) 
*                   Alen Stojanov      (astojanov@inf.ethz.ch)
*                   Gagandeep Singh    (gsingh@inf.ethz.ch)
*                   Markus Pueschel    (pueschel@inf.ethz.ch)
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see http://www.gnu.org/licenses/.
*/
//#include "stdafx.h"

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "tsc_x86.h"

using namespace std;

#define CYCLES_REQUIRED 1e7
#define REP 20
#define EPS (1e-3)
#define FREQ 2.7
#define TOLERANCE 1e-8
#define MAX_SIZE  300

/* prototype of the function you need to optimize */
typedef void(*comp_func)( const double*, int, double* );

//headers
double get_perf_score(comp_func f);
void register_functions();
double perf_test(comp_func f, int n);

//You can delcare your functions here
extern "C" void mil2_baseline(const double *hr_sphere_region, int n, double *directions_vectors_mil);
extern "C" void mil2_o1(const double *hr_sphere_region, int n, double *directions_vectors_mil);
//extern "C" void mil2_o2(const int *hr_sphere_region, int n, double *directions_vectors_mil);
extern "C" void dummy1(const double *hr_sphere_region, int n, double *directions_vectors_mil);
//extern "C" void dummy2(const int *hr_sphere_region, int n, double *directions_vectors_mil);
//extern "C" void dummy3(const int *hr_sphere_region, int n, double *directions_vectors_mil);

void add_function(comp_func f, const string& name, double flop);

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<double> funcFlops;
int numFuncs = 0;

template <typename T>
void rands(T * m, size_t n)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (size_t i = 0; i < n; ++i)
        m[i] = dist(gen) > 0.0 ? 1 : 0;
}

void build(int **a, int n)
{
    *a = static_cast<int *>(aligned_alloc(32, n * sizeof(int)));
    rands(*a, n);
}

void build(double **a, int n)
{
    *a = static_cast<double *>(aligned_alloc(32, n * sizeof(double)));
    rands(*a, n);
}

void destroy(void * m)
{
    free(m);
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    add_function(&mil2_baseline, "Base line", 3.25);
    add_function(&mil2_o1, "Base opt1", 3.25);
//    add_function(&mil2_o2, "Base opt2", 6.5);
//    add_function(&dummy1, "Dummy 1", 3.25);
//    add_function(&dummy2,   "Dummy 2", 3.25);
//    add_function(&dummy3,   "Dummy 3", 3.25);
}

bool checksum(const double* a, const double* b, int n) {

    for(int i = 0; i < n; i++) {
        if ( (b[i] < a[i] - TOLERANCE) || (b[i] > a[i] + TOLERANCE)) {
            return true;
        }
    }

    return false;
}

/*
* Main driver routine - calls register_funcs to get student functions, then
* tests all functions registered, and reports the best performance
*/
int main(int argc, char **argv)
{
    cout << endl << "Starting program. " << endl;
    double cycles;
    int i;

    register_functions();

    if (numFuncs == 0)
    {
        cout << endl;
        cout << "No functions registered - nothing for driver to do" << endl;
        cout << "Register functions by calling register_func(f, name)" << endl;
        cout << "in register_funcs()" << endl;

        return 0;
    }
    cout << numFuncs << " functions registered." << endl;

    //Check validity of functions. 
//    int n = 30;
    double *region;
    double* output;
    double* outputBaseline;

    for (int n = 20; n <= MAX_SIZE; n += 20) {
        cout << endl << "Testing size " << n << endl;

        // Compute with base line first.
        if (numFuncs > 1) {
            build(&region, n*n*n);
            build(&outputBaseline, 13);
            build(&output, 13);

            // First compute with baseline
            userFuncs[0](region, n, outputBaseline);

            for (i = 1; i < numFuncs; i++)
            {
                comp_func f = userFuncs[i];
                f(region, n, output);
                bool error = checksum(outputBaseline, output, 13);
                if (error)
                    cout << "ERROR: the results for function " << i << " are incorrect." << std::endl;
            }

            destroy(region);
            destroy(output);
            destroy(outputBaseline);
        }


        for (i = 0; i < numFuncs; i++)
        {
            cycles = perf_test(userFuncs[i], n);
            cout << endl << "** Running: " << funcNames[i] << " **" << endl;
            cout << "Runtime: " << cycles << " cycles" << endl;
            cout << "Performance: " << (1.0 * funcFlops[i] * n * n * n) / cycles << " flops per cycle" << endl;
//            cout << (1.0 * funcFlops[i] * n * n * n) / cycles << endl;
        }
    }

    return 0;
}


/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, const string& name, double flops)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flops);

    numFuncs++;
}

/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
double perf_test(comp_func f, int n)
{
    double cycles;
    long num_runs = 10;
    double multiplier = 1;
    myInt64 start, end;

    double *region;
    double* output;
    build(&region, n*n*n);
    build(&output, 13);


    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f(region, n, output);
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(region, n, output);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;

        cyclesList.push_back(cycles);
    }

    destroy(region);
    destroy(output);
    cyclesList.sort();
    cycles = cyclesList.front();
    return  cycles;
}


