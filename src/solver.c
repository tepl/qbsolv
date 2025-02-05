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
// Flips a given bit in the solution, and calculates the new energy.
//
// All the auxiliary information (flip_cost) is updated.
//
// @param old_energy The current objective function value
// @param bit is the bit to be flipped
// @param[in,out] solution inputs a current solution, flips the given bit
// @param qubo_size is the number of variables in the QUBO matrix
// @param qubo the QUBO matrix being solved
// @param[out] flip_cost The change in energy from flipping a bit
// @returns New energy of the modified solution
double evaluate_1bit(double old_energy, uint bit, int8_t *solution, uint qubo_size,
    double **qubo, double *flip_cost)
{
    double result = old_energy + flip_cost[bit];

    // Flip the bit and reverse its flip_cost
    solution[bit] = 1 - solution[bit];
    flip_cost[bit] = -flip_cost[bit];

    // Update the flip cost for all of the adjacent variables
    if (solution[bit] == 0 ) {
        // for rows ii up to bit, the flip_cost[ii] changes by qubo[bit][ii]
        // for columns ii from bit+1 and higher, flip_cost[ii] changes by qubo[ii][bit]
        // the sign of the change (positive or negative) depends on both solution[bit] and solution[ii].
        for (uint ii = 0; ii < bit; ii++)
            flip_cost[ii] += qubo[ii][bit]*(solution[ii]-!solution[ii]);

        for (uint ii = bit + 1; ii < qubo_size; ii++)
            flip_cost[ii] += qubo[bit][ii]*(solution[ii]-!solution[ii]);

    } else {
        // if solution[bit] was a 1 before, flip_cost[ii] changes in the other direction.
        for (uint ii = 0; ii < bit; ii++)
            flip_cost[ii] -= qubo[ii][bit]*(solution[ii]-!solution[ii]);

        for (uint ii = bit + 1; ii < qubo_size; ii++)
            flip_cost[ii] -= qubo[bit][ii]*(solution[ii]-!solution[ii]);

    }
    // proposed code to clean up some numerical noise
    //for ( uint ii=0; ii< qubo_size; ii++) 
    //      flip_cost[ii] = roundit(flip_cost[ii],10) ;   


    return result;
}


// Tries to improve the current solution Q by flipping single bits.
// It flips a bit whenever a bit flip improves the objective function value,
// terminating when a local optimum is found.
// It returns the objective function value for the new solution.
//
// This routine does not perform a full evalution of the the state or auxiliary
// information, it assumes it is already up to date.
//
// @param energy The current objective function value
// @param[in,out] solution inputs a current solution, modified by local search
// @param[in] size is the number of variables in the QUBO matrix
// @param[in] qubo the QUBO matrix being solved
// @param[out] flip_cost The change in energy from flipping a bit
// @param[in,out] bit_flips is the number of candidate bit flips performed in the entire algorithm so far
// @returns New energy of the modified solution
double local_search_1bit(double energy, int8_t *solution, uint qubo_size,
    double **qubo, double *flip_cost, int64_t *bit_flips)
{
    int kkstr = 0, kkend = qubo_size, kkinc;
    int *index;
    if (GETMEM(index, int, qubo_size) == NULL) BADMALLOC

    for (uint kk = 0; kk < qubo_size; kk++) {
        index[kk] = kk;
    }

    // The local search terminates at the local optima, so the moment we can't
    // improve with a single bit flip
    bool improve = true;
    while (improve) {
        improve = false;

        if (kkstr == 0) { // sweep top to bottom
            shuffle_index(index, qubo_size);
            kkstr = qubo_size - 1; kkinc = -1; kkend = 0;
        } else { // sweep bottom to top
            kkstr = 0; kkinc = 1; kkend = qubo_size; // got thru it backwards then reshuffle
        }

        for (int kk = kkstr; kk != kkend; kk = kk + kkinc) {
            uint bit = index[kk];
            (*bit_flips)++;
            if (energy + flip_cost[bit] > energy) {
                energy  = evaluate_1bit(energy, bit, solution, qubo_size, qubo, flip_cost);
                improve = true;
            }
        }
    }
    free(index);
    return energy;
}

// Performs a local Max search improving the solution and returning the last evaluated value
//
// Mostly the same as local_search_1bit, except it first evaluates the
// current solution and updates the auxiliary information (flip_cost)
//
// @param[in,out] solution inputs a current solution, modified by local search
// @param size is the number of variables in the QUBO matrix
// @param[out] flip_cost The change in energy from flipping a bit
// @param bit_flips is the number of candidate bit flips performed in the entire algorithm so far
// @returns New energy of the modified solution
double local_search(int8_t *solution, int qubo_size, double **qubo,
    double *flip_cost, int64_t *bit_flips)
{
    double energy;

    // initial evaluate needed before evaluate_1bit can be used
    energy = evaluate(solution, qubo_size, qubo, flip_cost);
    energy = local_search_1bit(energy, solution, qubo_size, qubo, flip_cost, bit_flips); // local search to polish the change
    return energy;
}

// This function is called by solve to execute a tabu search, This is THE Tabu search
//
// A tabu optimization algorithm tries to find an approximately maximal solution
// to a QUBO problem. Tabu is an optimization algorithm which performs single
// bit flips on a current solution in an attempt to improve it. For each
// candidate bit flip, the change in objective function value is examined.
// If this results in a new best solution, that change is accepted.
// If no bit flip results in a new best solution, we choose the best among
// the candidate bit flips.
//
// In order to avoid getting stuck in local optima, a list of "tabu" (not-allowed)
// moves is maintained (the vector "tabuK") After a bit has been flipped,
// it cannot be flipped again for another "nTabu" moves. The algorithm terminates
// after sufficiently many bit flips without improvment.
//
// @param[in,out] solution inputs a current solution and returns the best solution found
// @param[out] best stores the best solution found during the algorithm
// @param qubo_size is the number of variables in the QUBO matrix
// @param qubo is the QUBO matrix to be solved
// @param flip_cost is the impact vector (the change in objective function value that results from flipping each bit)
// @param bit_flips is the number of candidate bit flips performed in the entire algorithm so far
// @param iter_max is the maximum size of bit_flips allowed before terminating
// @param TabuK stores the list of tabu moves
// @param target Halt if this energy is reached and TargetSet is true
// @param target_set Do we have a target energy at which to terminate
// @param index is the order in which to perform candidate bit flips (determined by flip_cost).
double tabu_search(int8_t *solution, int8_t *best, uint qubo_size, double **qubo,
    double *flip_cost, int64_t *bit_flips, int64_t iter_max,
    int *TabuK, double target, bool target_set, int *index, int nTabu)
{

    uint      last_bit = 0; // Track what the previously flipped bit was
    bool      brk; // flag to mark a break and not a fall-thru of the loop
    double    best_energy; // best solution so far
    double    Vlastchange; // working solution variable
    double    sign;
    int64_t thisIter;
    int64_t increaseIter;
    int       numIncrease = 900;
    double    howFar;

    // setup nTabu
    // these nTabu numbers might need to be adjusted to work correctly
    if (nTabu == 0 ) {   // nTabu not specified on call
        if (Tlist_ != -1) {
            nTabu = MIN(Tlist_, (int)qubo_size + 1 ); // tabu use set tenure
        } else {
            if      (qubo_size <  20)  nTabu = 5;
            else if (qubo_size < 100)  nTabu = 10;
            else if (qubo_size < 250)  nTabu = 12;
            else if (qubo_size < 500)  nTabu = 13;
            else if (qubo_size < 1000) nTabu = 21;
            else if (qubo_size < 2500) nTabu = 29;
            else if (qubo_size < 8000) nTabu = 34;
            else /*qubo_size >= 8000*/ nTabu = 35;
        }
    }

    if ( findMax_ ) {
        sign = 1.0;
    } else {
        sign = -1.0;
    }

    best_energy  = local_search(solution, qubo_size, qubo, flip_cost, bit_flips);
    val_index_sort(index, flip_cost, qubo_size); // Create index array of sorted values
    thisIter     = iter_max - (*bit_flips);
    increaseIter = thisIter / 2;
    Vlastchange  = best_energy;

    for (uint i = 0; i < qubo_size; i++) best[i] = solution[i]; // copy the best solution so far
    for (uint i = 0; i < qubo_size; i++) TabuK[i] = 0; // zero out the Tabu vector

    int kk, kkstr = 0, kkend = qubo_size, kkinc;
    while (*bit_flips < iter_max) {
        // best solution in neighbour, initialized most negative number
        double neighbour_best = BIGNEGFP;
        brk = false;
        if ( kkstr == 0 ) { // sweep top to bottom
            kkstr = qubo_size - 1; kkinc = -1; kkend = 0;
        } else { // sweep bottom to top
            kkstr = 0; kkinc = 1; kkend = qubo_size;
        }

        for (kk = kkstr; kk != kkend; kk = kk + kkinc) {
            uint bit = index[kk];
            if (TabuK[bit] != (int8_t)0 ) continue;
            {
                (*bit_flips)++;
                double new_energy = Vlastchange + flip_cost[bit]; //  value if Q[k] bit is flipped
                if (new_energy > best_energy) {
                    brk         = true;
                    last_bit    = bit;
                    new_energy  = evaluate_1bit(Vlastchange, bit, solution, qubo_size, qubo, flip_cost); // flip the bit and fix tables
                    Vlastchange = local_search_1bit(new_energy, solution, qubo_size, qubo, flip_cost, bit_flips); // local search to polish the change
                    val_index_sort_ns(index, flip_cost, qubo_size); // update index array of sorted values, don't shuffle index
                    best_energy = Vlastchange;

                    for (uint i = 0; i < qubo_size; i++) best[i] = solution[i]; // copy the best solution so far

                    if (target_set) {
                        if (Vlastchange >= (sign * target)) {
                            break;
                        }
                    }
                    howFar = ((double)(iter_max - (*bit_flips)) / (double)thisIter);
                    if (Verbose_ > 3) {
                        printf("Tabu new best %lf ,K=%d,iteration = %"LONGFORMAT", %lf, %d\n",
                            best_energy * sign, last_bit,(int64_t) (*bit_flips), howFar, brk );
                    }
                    if ( howFar < 0.80  && numIncrease > 0 ) {
                        if (Verbose_ > 3) {
                            printf("Increase Itermax %"LONGFORMAT", %"LONGFORMAT"\n",  iter_max, 
                                     (iter_max + increaseIter));
                        }
                        iter_max  += increaseIter;
                        thisIter += increaseIter;
                        numIncrease--;
                    }
                    break;
                }
                // Q vector unchanged
                if (new_energy > neighbour_best) { // check for improved neighbour solution
                    last_bit = bit;   // record position
                    neighbour_best = new_energy;   // record neighbour solution value
                }
            }
        }

        if (target_set) {
            if (Vlastchange >= (sign * target)) {
                break;
            }
        }
        if ( !brk ) { // this is the fall-thru case and we haven't tripped interior If V> VS test so flip Q[K]
            Vlastchange = evaluate_1bit(Vlastchange, last_bit, solution, qubo_size, qubo, flip_cost);
        }

        uint i;
        for (i = 0; i < qubo_size; i++) TabuK [i] =  MAX(0, TabuK[i] - 1);

        // add some asymmetry
        if (solution[qubo_size-1] == 0) {
            TabuK[last_bit] = nTabu + 1;
        } else {
            TabuK[last_bit] = nTabu - 1;
        }
    }

    // copy over the best solution
    for (uint i = 0; i < qubo_size; i++) solution[i] = best[i];

    // ok, we are leaving Tabu, we can do a for-sure clean-up run of evaluate, to be sure we
    // return the true evaluation of the function (given that we only do this a handful of times)
    double final_energy = evaluate(solution, qubo_size, qubo, flip_cost);

    // Create index array of sorted values
    val_index_sort(index, flip_cost, qubo_size);
    return final_energy;
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

        // Go over all other (non-extracted) variables
        // from the highest until we reach the current variable
        for (int j = qubo_size - 1; j > variable; j--) {
            if ( j == Icompress[ji] ) {
                // Found a sub_qubo element, skip it, watch for the next one
                ji--;
            } else {
                clamp += qubo[variable][j] * solution[j];
            }
        }

        // Go over all other (non-extracted) variables
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

// solv_submatrix() performs QUBO optimization on a subregion.
// In this function the subregion is optimized using tabu_search() rather than using the D-Wave hardware.
//
// @param[in,out] solution inputs a current solution and returns the best solution found
// @param[out] best stores the best solution found during the algorithm
// @param qubo_size is the number of variables in the QUBO matrix
// @param qubo is the QUBO matrix to be solved
// @param flip_cost is the impact vector (the change in objective function value that results from flipping each bit)
// @param bit_flips is the number of candidate bit flips performed in the entire algorithm so far
// @param TabuK stores the list of tabu moves
// @param index is the order in which to perform candidate bit flips (determined by Qval).
double solv_submatrix(int8_t *solution, int8_t *best, uint qubo_size, double **qubo,
    double *flip_cost, int64_t *bit_flips, int *TabuK, int *index)
{
    int nTabu;
    int64_t iter_max = (*bit_flips) + (int64_t)MAX((int64_t)3000, (int64_t)20000 * (int64_t)qubo_size);
    if      (qubo_size <  20)  nTabu = 5;
    else if (qubo_size < 100)  nTabu = 10;
    else if (qubo_size < 250)  nTabu = 12;                                                                
    else if (qubo_size < 500)  nTabu = 13;                                                                
    else if (qubo_size < 1000) nTabu = 21;
    else if (qubo_size < 2500) nTabu = 29;                                                                
    else if (qubo_size < 8000) nTabu = 34;                                                                
    else /*qubo_size >= 8000*/ nTabu = 35;
    
    return tabu_search(solution, best, qubo_size, qubo, flip_cost,
        bit_flips, iter_max, TabuK, Target_, false, index, nTabu);
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
    //get scratch memory needed by tabu solver
    int8_t sub_solution[subMatrix], Best[subMatrix];
    double **sub_qubo, flip_cost[subMatrix];
    int64_t bit_flips;
    int *TabuK, *index;

    sub_qubo = (double**)malloc2D(qubo_size, qubo_size, sizeof(double));
    if (GETMEM(index, int, subMatrix) == NULL) BADMALLOC
    if (GETMEM(TabuK, int, subMatrix) == NULL) BADMALLOC

    reduce(Icompress, qubo, subMatrix, qubo_size, sub_qubo, solution, sub_solution);
    // solve
    if (Verbose_ > 3) { 
        printf("\nBits set before solver "); 
        for (int j = 0; j < subMatrix; j++) printf("%d", solution[Icompress[j]]); 
    }
    if ( UseDwave_ ) {
        dw_solver(sub_qubo, subMatrix, sub_solution);
        int64_t sub_bit_flips=0;  //  run a local search with higher precision than the Dwave
        local_search(sub_solution, subMatrix, sub_qubo,flip_cost, &sub_bit_flips);
    }else {
        bit_flips=0;
        for (int i = 0;i<subMatrix;i++ ) {
            TabuK[i]=0;
            index[i]=i;
            sub_solution[i]=solution[Icompress[i]];
            Best[i]=solution[Icompress[i]];
        }
        solv_submatrix(sub_solution, Best, subMatrix, sub_qubo, flip_cost, &bit_flips, TabuK, index);
    } 
                    //char subqubofile[sizeof "subqubo10000.qubo"]; // modification to write out subqubos
                    //sprintf(subqubofile,"subqubo%05ld.qubo",numPartCalls);
                    //write_qubo(sub_qubo,subMatrix,subqubofile);
    // projection
    if (Verbose_ > 3) { 
        printf("\nBits set after solver  "); 
        for (int j = 0; j < subMatrix; j++) printf("%d", sub_solution[j]); 
        printf("\n");
    }
    for (int j = 0; j < subMatrix; j++) {
        int bit = Icompress[j];
        if (solution[bit] != sub_solution[j] ) change++;
        solution[bit] = sub_solution[j];
    }
    free( sub_qubo );
    free( TabuK );
    free( index );
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
// When none of the variables in any of the subregions change, a new solution
// is chosen based on randomizing those variables.
//
// After Pchk = 8 iterations with no improvement, the algorithm is
//   completely restarted with a new random solution.
// After nRepeats iterations with no improvement, the algorithm terminates.
//
// @param qubo The QUBO matrix to be solved
// @param qubo_size is the number of variables in the QUBO matrix
// @param nRepeats is the number of iterations without improvement before giving up
void solve(double **qubo, const int qubo_size, int nRepeats, int8_t **solution_list, double *energy_list, int *solution_counts, int *Qindex, int QLEN)
{
    double energy;
    double best_energy;
    int    start_;
    long   numPartCalls = 0;
    double sign = findMax_ ? 1.0 : -1.0;
    double **mat;
    int8_t *solution;
    int8_t *Qbest;
    int8_t  *matSol;

    start_ = clock();

    // allocate solution vector
    if (GETMEM(solution, int8_t, qubo_size) == NULL) BADMALLOC
    for (int i = 0; i < qubo_size; i++) {
        solution[i]=0;
    }

    // allocate solutions storage
    int num_nq_solutions=0;
    for (int i = 0; i < QLEN+1 ; i++) {
        energy_list[i]     = BIGNEGFP;
        solution_counts[i] = 0;
        for (int j = 0; j < qubo_size; j++ ) {
            solution_list[i][j] = 0;
        }
    }

    // exit if qubo_size is large then SubMatrix_
    if (qubo_size > SubMatrix_) {
        printf("Qubo_size > SubMatrix_, exiting...\n");
        return;

    // or start solving
    } else {

        // allocate temporary arrays
        mat = (double**)malloc2D(SubMatrix_, SubMatrix_, sizeof(double) );
        if (GETMEM(matSol, int8_t, SubMatrix_) == NULL) BADMALLOC

        // loop over repeats
        for (int i = 0; i < nRepeats; i++) {

            // fill out matrix
            for (int i = 0; i < SubMatrix_; i++) {
                for (int j = 0; j < SubMatrix_; j++)
                    if( i < qubo_size && j < qubo_size ) {
                        mat[i][j] = qubo[i][j];
                    } else {
                        mat[i][j] = 0.0;
                    }
            }

            // run D-Wave solver
            dw_solver(mat, SubMatrix_, matSol);

            // copy solution
            for (int i = 0; i < qubo_size; i++) {
              solution[i] = matSol[i];
            }

            // evaluate energy
            energy = Simple_evaluate(solution, qubo_size, qubo);

            // store solution
            manage_solutions(solution, solution_list, energy, energy_list, solution_counts, Qindex, QLEN, qubo_size, &num_nq_solutions);
        }
    }

    // print best solution
    Qbest = &solution_list[Qindex[0]][0];
    best_energy = energy_list[Qindex[0]];
    print_output(qubo_size, Qbest, numPartCalls, best_energy * sign, CPSECONDS);

    free(mat);
    free(matSol);
    free(solution);

    return;
}
