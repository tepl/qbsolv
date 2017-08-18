/*
 *
 Copyright 2017 D-Wave Systems Inc

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

#include "include.h"
#include "extern.h"
#include <math.h>

// This function Simply evaluates the objective function for a given solution.
//  call this to double check that a solution = energy
// @param solution a current solution
// @param qubo_size the number of variables in the QUBO matrix
// @param qubo the QUBO matrix being solved
// @returns Energy of solution evaluated by qubo
double Simple_evaluate(int8_t *solution, uint qubo_size, double **qubo)
{
    double result = 0.0;

    for (uint ii = 0; ii < qubo_size; ii++) {
        double row_sum = 0.0;
        double col_sum = 0.0;

        // qubo an upper triangular matrix, so start right of the diagonal
        // for the rows, and stop at the diagonal for the columns
        for (uint jj = ii + 1; jj < qubo_size; jj++)
            row_sum += qubo[ii][jj] * (double)solution[jj];

        for (uint jj = 0; jj < ii; jj++)
            col_sum += qubo[jj][ii] * (double)solution[jj];

        if (solution[ii] == 1) result += row_sum + qubo[ii][ii];
    }

    return result;
}
// This function evaluates the objective function for a given solution.
//
// It is called when the search is starting over, such as after a projection
// in outer loop of solver.
//
// @param solution a current solution
// @param qubo_size the number of variables in the QUBO matrix
// @param qubo the QUBO matrix being solved
// @param[out] flip_cost The change in energy from flipping a bit
// @returns Energy of solution evaluated by qubo
double evaluate(int8_t *solution, uint qubo_size, double **qubo, double *flip_cost)
{
    double result = 0.0;

    for (uint ii = 0; ii < qubo_size; ii++) {
        double row_sum = 0.0;
        double col_sum = 0.0;

        // qubo an upper triangular matrix, so start right of the diagonal
        // for the rows, and stop at the diagonal for the columns
        for (uint jj = ii + 1; jj < qubo_size; jj++)
            row_sum += qubo[ii][jj] * (double)solution[jj];

        for (uint jj = 0; jj < ii; jj++)
            col_sum += qubo[jj][ii] * (double)solution[jj];

        // If the variable is currently 1, then by flipping it we lose
        // what it is currently contributing (so we negate the contribution),
        // when it is currently false, gain that ammount by flipping
        if (solution[ii] == 1) {
            result += row_sum + qubo[ii][ii];
            flip_cost[ii] = -(row_sum + col_sum + qubo[ii][ii]);
        } else {
            flip_cost[ii] =  (row_sum + col_sum + qubo[ii][ii]);
        }
    }

    return result;
}

// reduce() computes a subQUBO (val_s) from large QUBO (val)
// for the purposes of optimizing on a subregion (a subset of the variables).
// It does this by fixing all variables outside the subregion to their current values,
// and adding the influence of the fixed variables to the linear (diagonal) terms in the subQUBO.
//
// @param Icompress is the list of variables in the subregion that will be extracted
// @param qubo is the large QUBO matrix to be solved
// @param sub_qubo_size is the number of variable in the subregion
// @param qubo_size is the number of variables in the large QUBO matrix
// @param[out] sub_qubo is the returned subQUBO
// @param[out] sub_solution is a current solution on the subQUBO
void reduce(int *Icompress, double **qubo, uint sub_qubo_size, uint qubo_size,
    double **sub_qubo, int8_t *solution, int8_t *sub_solution)
{
    // clean out the subMatrix
    for (uint i = 0; i < sub_qubo_size; i++) { // for each column
        for (uint j = 0; j < sub_qubo_size; j++)
            sub_qubo[i][j] = 0.0; // for each row
    }

    // fill the subMatrix
    for (uint i = 0; i < sub_qubo_size; i++) { // for each column
        sub_solution[i] = solution[Icompress[i]];
        for (uint j = i; j < sub_qubo_size; j++) { // copy row
            sub_qubo[i][j] = qubo[Icompress[i]][Icompress[j]];
        }
    }

    // The remainder of the function is clamping the sub_qubo to the
    // solution state surrounding it.
    // Go over every variable that we are extracting
    for (uint sub_variable = 0; sub_variable < sub_qubo_size; sub_variable++) {
        // Get the global index of the current variable
        int variable = Icompress[sub_variable];
        double clamp = 0;

        // this will keep track of the index of the next sub_qubo component,
        // we don't include those in the clamping
        int ji = sub_qubo_size - 1;

        // Go over all other (non-extracted) variable
        // from the highest until we reach the current variable
        for (int j = qubo_size - 1; j > variable; j--) {
            if ( j == Icompress[ji] ) {
                // Found a sub_qubo element, skip it, watch for the next one
                ji--;
            } else {
                clamp += qubo[variable][j] * solution[j];
            }
        }

        // Go over all other (non-extracted) variable
        // from zero until we reach the current variable
        ji = 0;
        for (int j = 0; j < variable + 1; j++) {
            if ( j == Icompress[ji] ) {
                // Found a sub_qubo element, skip it, watch for the next one
                ji++;
            } else {
                clamp += qubo[j][variable] * solution[j];
            }
        }

        // Now that we know what the effects of the non-extracted variables
        // are on the sub_qubo we include it by adding it as a linear
        // bias in the sub_qubo (a diagonal matrix entry)
        sub_qubo[sub_variable][sub_variable] += clamp ;
    }
    return;
}

// reduce_solv_projection reduces from a submatrix solves the QUBO projects the solution and 
//      returns the number of changes
// @param Icompress index vector , ordered lowest to highest, of the row/columns to extract subQubo
// @param qubo is the QUBO matrix to extract from
// @param qubo_size is the number of variables in the QUBO matrix                                                                    
// @param subMatrix is the size of the subMatrix to create and solve 
// @param @param[in,out] solution inputs a current solution and returns the projected solution 
// @param[out] stores the new, projected solution found during the algorithm

int reduce_solve_projection( int *Icompress, double **qubo, int qubo_size, int subMatrix, int8_t *solution )
{
    int change=0;
    int8_t sub_solution[subMatrix];
    double **sub_qubo;

    sub_qubo = (double**)malloc2D(qubo_size, qubo_size, sizeof(double));

    reduce(Icompress, qubo, subMatrix, qubo_size, sub_qubo, solution, sub_solution);

    // solve
    dw_solver(sub_qubo, subMatrix, sub_solution);

    // projection
    for (int j = 0; j < subMatrix; j++) {
        int bit = Icompress[j];
        if (solution[bit] != sub_solution[j] ) change++;
        solution[bit] = sub_solution[j];
    }
    return change;
}

// Entry into the overall solver from the main program
//
// It is the main function for solving a quadratic boolean optimization problem.
//
// The algorithm alternates between:
//   1) performing a global tabu search (tabu_search()) from the current solution
//   2) selecting subregions and optimizing on each of those subregions with
//      all other variables clamped (solv_submatrix()).
//
// Subregions are chosen based on increasing impact ("flip_cost"). Impact
// is defined as the change in objective function value that occurs by flipping a bit.
// The first subregion is chosen as the N least impactful variables, and so on.
// When none of the variables any in the subregions change, a new solution
// is chosen based on randomizing those variables.
//
// After Pchk = 8 iterations with no improvement, the algorithm is
//   completely restarted with a new random solution.
// After nRepeats iterations with no improvement, the algorithm terminates.
//
// @param qubo The QUBO matrix to be solved
// @param qubo_size is the number of variables in the QUBO matrix
// @param nRepeats is the number of iterations without improvement before giving up
void solve(double **qubo, const int qubo_size, int nRepeats)
{
    double    *flip_cost, energy;
    int       *index, start_;
    int8_t    *solution;
    int       *Icompress;
    int        subMatrix    = SubMatrix_;
    long       numPartCalls = 0;
    double    *energies;
    int8_t    *best_solution;
    double     best_energy;
    int        best_repeat;

    start_ = clock();

    // Allocate memory
    if (GETMEM(solution, int8_t, qubo_size) == NULL) BADMALLOC
    if (GETMEM(best_solution, int8_t, qubo_size) == NULL) BADMALLOC
    if (GETMEM(flip_cost, double, qubo_size) == NULL) BADMALLOC
    if (GETMEM(index, int, qubo_size) == NULL) BADMALLOC
    if (GETMEM(Icompress, int, qubo_size) == NULL) BADMALLOC
    if (GETMEM(energies, double, nRepeats) == NULL) BADMALLOC

    // QUBO fits into one subQUBO
    if ( qubo_size <= subMatrix ) {

        dw_solver(qubo, qubo_size, best_solution);

        best_energy = Simple_evaluate(best_solution,  qubo_size, qubo);

    }

    // QUBO does not fit into one subQUBO
    else {

        // Setup l_max
        int  l_max = qubo_size - SubMatrix_;

        // Setup best energy
        best_energy = BIGNEGFP;

        // Use random solution to initialize flip_cost
        randomize_solution(solution, qubo_size);
        energy = evaluate(solution, qubo_size, qubo, flip_cost);

        // Start outer loop
        for (int RepeatPass = 0; RepeatPass < nRepeats; RepeatPass++){

            // Sort by impact
            val_index_sort(index, flip_cost, qubo_size);

            // Submatrix passes
            int change=0;
            for (int l = 0; l < l_max; l += subMatrix) {
                if (Verbose_ > 3) printf("Submatrix starting at backbone %d\n", l);

                for (int i = l, j = 0; i < l + subMatrix; i++) {
                    Icompress[j++] = index[i]; // create compression index
                }
                index_sort(Icompress, subMatrix, true); // sort it for effective reduction

                change=change+reduce_solve_projection( Icompress, qubo, qubo_size, subMatrix, solution );
                numPartCalls++;
            }

            // Compute energy and flip_cost
            energy = evaluate(solution, qubo_size, qubo, flip_cost);

            // Save energy
            energies[RepeatPass] = energy;

            // Update best solution
            if(energy > best_energy){
                best_energy = energy;
                best_repeat = RepeatPass;
                for (int i = 0; i < qubo_size; i++) {
                    best_solution[i] = solution[i];
                }
            }

        } // End outer loop

    }

    // Print results
    double sign = findMax_ ? 1.0 : -1.0;
    print_output(qubo_size, best_solution, numPartCalls, sign * best_energy, CPSECONDS);

    // Print passes
    if( qubo_size > subMatrix ) {
        for (int RepeatPass = 0; RepeatPass < nRepeats; RepeatPass++){
            fprintf(outFile_, "RepeatPass %d %f\n", RepeatPass, sign * energies[RepeatPass]);
        }
        fprintf(outFile_, "Best %d %f\n", best_repeat, sign * energies[best_repeat]);
    }

    // Write matrix if needed
    if (WriteMatrix_)
        print_solution_and_qubo(solution, qubo_size, qubo);

    // Free memory
    free(solution);
    free(best_solution);
    free(flip_cost);
    free(index);
    free(Icompress);
    free(energies);
    free(qubo); 
    return;

}
