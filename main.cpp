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
#include <stdio.h>
#include <reader.h>
#include "tsc_x86.h"
#include "kernel.h"

using namespace std;

#define CYCLES_REQUIRED 1e7
#define REP 3
#define EPS (1e-3)
#define CLOCK_FREQ 2.7

/* prototype of the function you need to optimize */
typedef void(*comp_func)(int*, int*, int*, double*, double*, double *, int);

//headers
double get_perf_score(comp_func f);
void register_functions();
double perf_test(comp_func f, string desc, long flops, int);

//You can delcare your functions here
//void kernel_basic(int* sphere, double* ptrHighRes, double* ptrLowRes, double* rotation_matrix, double* ptrEvecOut);

void add_function(comp_func f, string name, long flop);

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<long> funcFlops;
int numFuncs = 0;

void build(int** sphere, int** ptrHighRes, int** ptrLowRes, double** rotation_matrix, double** ptrEvecOut, double **ptrEvalsOut, int index)
{
    init(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, index);
}

void destroy(int* sphere, int* ptrHighRes, int* ptrLowRes, double* rotation_matrix, double* ptrEvecOut, double *ptrEvalsOut, bool generate_ground_truth)
{
    deInit(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, generate_ground_truth);
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions()
{
    // TODO: Add correct number of flops
    add_function(&kernel_basic, "Base kernel", 1572181 + 126347988814);
    // Add your functions here. Don't modify the number of flops parameter.
}


/*
 * test evals and evecs against ground truth file
 *
 */
double checksum(double (*evals)[3], double (*evecs)[3][3], int n) 
{
    FILE *fd = fopen("ground_truth.txt", "r");

    double evecs_true[3][3], evals_true[3];
    int index[3];
    double error = 0;

    for (int i=0; i<n; i++) {
        /* read in ground truth */
        fscanf(fd, 
            "%d, %d, %d, " /* index into LR image */
            "%lf, %lf, %lf, "   /* e vals */
            "%lf, %lf, %lf, "   /* e vecs*/
            "%lf, %lf, %lf, "
            "%lf, %lf, %lf\n",
            &index[0], &index[1], &index[2],
            &evals_true[0], &evals_true[1], &evals_true[2],
            &evecs_true[0][0], &evecs_true[0][1], &evecs_true[0][2], 
            &evecs_true[1][0], &evecs_true[1][1], &evecs_true[1][2],
            &evecs_true[2][0], &evecs_true[2][1], &evecs_true[2][2]);
        
        /* assumption: index already matches */
        
        /* accumulate errors */
        for (int j=0; j<3; j++) {
            error += fabs(evals_true[j] - evals[i][j]);
            for (int k=0; k<3; k++) {
                error += fabs(evecs_true[k][j] - evecs[i][k][j]);
            }
        }
    }

    fclose(fd);
    return error;
}

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, string name, long flops)
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
double perf_test(comp_func f, string desc, long flops, int input)
{
    double cycles = 0.;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;

    int* sphere;
    int* ptrHighRes;
    int* ptrLowRes;
    double* rotation_matrix;
    double* ptrEvecOut;
    double* ptrEvalsOut;

    build(&sphere, &ptrHighRes, &ptrLowRes, &rotation_matrix, &ptrEvecOut, &ptrEvalsOut, input);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, input);
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
            f(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, input);
        }
        end = stop_tsc(start);

        cycles = ((double)end) / num_runs;

        cyclesList.push_back(cycles);
    }

    destroy(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, false);
    cyclesList.sort();
    cycles = cyclesList.front();
    return cycles;
}

/*
* Main driver routine - calls register_funcs to get student functions, then
* tests all functions registered, and reports the best performance
*/
int main(int argc, char **argv)
{
    cout << "Starting program. ";
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

    /* generate the ground truth data. assumption: userFuncs[0] generates ground truth */

    int*    sphere;
    int* ptrHighRes;
    int* ptrLowRes;
    double* rotation_matrix;
    double* ptrEvecOut;
    double* ptrEvalsOut;

    build(&sphere, &ptrHighRes, &ptrLowRes, &rotation_matrix, &ptrEvecOut, &ptrEvalsOut, 0);
    userFuncs[0](sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, 0);
    destroy(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, true);

    for (i =0; i < numFuncs; i++) {
        build(&sphere, &ptrHighRes, &ptrLowRes, &rotation_matrix, &ptrEvecOut, &ptrEvalsOut, 0);

        comp_func f = userFuncs[i];
        f(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, 0);
        double error = checksum((double (*)[3]) ptrEvalsOut, 
                                (double (*)[3][3]) ptrEvecOut, 
                                LOW_RES_SIZE(0));
        if (error > EPS) {
            cout << "ERROR: the results for function " << i << " are incorrect -> ";
            printf(" %g\n", error);
        } else {
            printf("userFuncs[%d] is correct!\n", i);
        }
        destroy(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, false);
    }

    for (i = 0; i < numFuncs; i++)
    {
        /* Iterate through all inputs */
        for (int j = 0; j < NUMBER_OF_INPUTS; ++j) {
            cout << endl << "*** Running " << funcNames[i] << " with input " << j << " ***" << endl;
            cycles = perf_test(userFuncs[i], funcNames[i], funcFlops[i], j);
            cout << "Runtime:     " << cycles << " cycles" << endl;
            cout << "Runtime:     " << cycles / CLOCK_FREQ / 1e6 << " msec" << endl;
            cout << "Performance: " << funcFlops[i] / cycles << " flops/cycle" << endl;
        }
    }
    cout << endl;

    free(ptrHighRes);
    return 0;
}

void simple_main () {
    int* sphere;
    int* ptrHighRes;
    int* ptrLowRes;
    double* rotation_matrix;
    double* ptrEvecOut;
    double* ptrEvalsOut;

    init(&sphere, &ptrHighRes, &ptrLowRes, &rotation_matrix, &ptrEvecOut, &ptrEvalsOut, 0);
    kernel_basic(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, 0);
    deInit(sphere, ptrHighRes, ptrLowRes, rotation_matrix, ptrEvecOut, ptrEvalsOut, false);

    printf("\nDone\n");
}
