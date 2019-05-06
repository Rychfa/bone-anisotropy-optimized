/**    ____                   __  ___
*     / __ )____  ____  ___  /  |/  /___ _____
*    / __  / __ \/ __ \/ _ \/ /|_/ / __ `/ __ \
*   / /_/ / /_/ / / / /  __/ /  / / /_/ / /_/ /
*  /_____/\____/_/ /_/\___/_/  /_/\__,_/ .___/
*                                     /_/
*  Bone Anisotropy Mapping.
*
*  Copyright (C) 2019  Jarunan Panyasantisuk (jarunan@ethz.ch)
*                      Rajan Gill            (rgill@ethz.ch)
*                      Ryan Cherifa          (rcherifa@student.ethz.ch)
*                      Joao Rivera           (hector.rivera@inf.ethz.ch)
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
*  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include "tsc_x86.h"
#include "kernel.h"

using namespace std;

#define CYCLES_REQUIRED 1e7
#define REP 25
#define EPS (1e-3)

/* prototype of the function you need to optimize */
typedef void(*comp_func)(int*, int*, int*, double*, double*);

//headers
double get_perf_score(comp_func f);
void register_functions();
double perf_test(comp_func f, string desc, int flops);

//You can delcare your functions here
//void kernel_basic(int* sphere, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut);

void add_function(comp_func f, string name, int flop);

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;

void build(int** sphere, int** ptrHighRes, int** ptrLowRes, double** rotation_matrix, double** ptrEvecOut)
{
    init(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);
}

void destroy(int* sphere, int* ptrHighRes, int* ptrLowRes, double* rotation_matrix, double* ptrEvecOut)
{
    deInit(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    // TODO: Add correct number of flops
    add_function(&kernel_basic, "Base kernel", 100);
    // Add your functions here. Don't modify the number of flops parameter.
}

// TODO: Add proper checksum
//double checksum(double *A, double *B, double *C, int n) {
//    double *w;
//    double *Cw;
//    double *Bw;
//    double *ABw;
//    build(&w, NR, 1);
//    build(&Cw, MR, 1);
//    build(&Bw, n, 1);
//    build(&ABw, MR, 1);
//
//    mvm(C, NR, 1, w, Cw, MR, NR);
//    mvm(B, NR, 1, w, Bw, n, NR);
//    mvm(A, 1, MR, Bw, ABw, MR, n);
//
//    double nrm_sqr = 0.0;
//    for(int i = 0; i < MR; i++) {
//        nrm_sqr += (Cw[i] - ABw[i]) * (Cw[i] - ABw[i]);
//    }
//
//    destroy(w);
//    destroy(Cw);
//    destroy(Bw);
//    destroy(ABw);
//
//    return nrm_sqr;
//}

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, string name, int flops)
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
double perf_test(comp_func f, string desc, int flops)
{
    double cycles = 0.;
    long num_runs = 5;
    double multiplier = 1;
    myInt64 start, end;

    int*    sphere;
    int* ptrHighRes;
    int* ptrLowRes;
    double* rotation_matrix;
    double* ptrEvecOut;

    build(&sphere, &ptrHighRes, &ptrLowRes, &rotation_matrix, &ptrEvecOut);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);
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
            f(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;

        cyclesList.push_back(cycles);
    }

    destroy(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);
    cyclesList.sort();
    cycles = cyclesList.front();
    return  cycles;
//    return  (1.0 * flops) / cycles;
}

/*
* Main driver routine - calls register_funcs to get student functions, then
* tests all functions registered, and reports the best performance
*/
int main(int argc, char **argv)
{
    cout << "Starting program. ";
    double perf;
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

    //
    // TODO: Implement validity of functions by comparing output with
    //  baseline before performance test.
    //

    for (i = 0; i < numFuncs; i++)
    {
        perf = perf_test(userFuncs[i], funcNames[i], funcFlops[i]);
        cout << "Running: " << funcNames[i] << endl;
        cout << perf << " cycles" << endl;
    }
    cout << endl;

    return 0;
}

int simple_main () {
    int* sphere;
    int* ptrHighRes;
    int* ptrLowRes;
    double* rotation_matrix;
    double* ptrEvecOut;

    init(&sphere, &ptrHighRes, &ptrLowRes, &rotation_matrix, &ptrEvecOut);
    kernel_basic(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);
    deInit(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut);

    printf("\nDone\n");

    return 0;
}
