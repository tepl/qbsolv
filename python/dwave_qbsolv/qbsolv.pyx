from libc.stdint cimport int8_t, int64_t
from libc.stdio cimport stdout
from libc.stdlib cimport malloc, free, srand
from cqbsolv cimport Verbose_, SubMatrix_, UseDwave_, algo_, outFile_, Time_, Tlist_, numsolOut_, WriteMatrix_, TargetSet_, findMax_
from cqbsolv cimport solve, malloc2D
import random
import sys

ENERGY_IMPACT = 0
SOLUTION_DIVERSITY = 1

PY2 = sys.version_info[0] == 2
if PY2:
    iteritems = lambda d: d.iteritems()
else:
    iteritems = lambda d: d.items()

def run_qbsolv(Q, repeats=50, 
               seed=17932241798878, verbosity=-1, algorithm=None, timeout=2592000):
    """Runs qbsolv.
    
    Args:
        Q (dict): A dictionary defining the QUBO. Should be of the form
            {(u, v): bias} where u, v are variables and bias is numeric.
        repeats (int, optional): Determines the number of times to
            repeat the main loop in qbsolv after determining a better
            sample. Default 50.
        seed (int, optional): Random seed. Default None.
        algorithm (int, optional): Which algorithm to use. Defaut
            None. Algorithms numbers can be imported from module
            under the names ENERGY_IMPACT and SOLUTION_DIVERSITY.

    Returns:
        (list, list, list): (samples, energies, counts) where samples is
        a list of dics, energies is the energy for each sample and counts
        is the number of times each sample was found by qbsolv.

    NB: relies on variables in Q being index valued and nonnegative, but not checked at this point
    """

    # qbsolv relies on a number of global variables that need to be set before the algorithm can
    # be run, so let's go ahead and set them here.

    # some of the global variables are changed based on input parameters
    global Verbose_
    Verbose_ = verbosity

    global algo_
    cdef int n_solutions = 20  # the maximimum number of solutions returned
    if algorithm is None or algorithm == ENERGY_IMPACT:
        algo_[0] = "o"
        algo_[1] = 0
        # n_solutions = 20
    elif algorithm == SOLUTION_DIVERSITY:
        algo_[0] = "d"
        algo_[1] = 0
        n_solutions = 70
    else:
        raise ValueError('unknown algorithm given')

    if not isinstance(timeout, int) or timeout <= 0:
        raise ValueError("'timeout' must be a positive integer")
    global Time_
    Time_ = timeout # the maximum runtime of the algorithm in seconds before timeout (2592000 = a months worth of seconds)

    # other global variables are fixed for our purposes here
    global UseDwave_
    UseDwave_ = False

    global outFile_
    outFile_ = stdout

    global SubMatrix_
    SubMatrix_ = 46

    global Tlist_
    Tlist_ = -1

    global numsolOut_
    numsolOut_ = 0

    global WriteMatrix_
    WriteMatrix_ = False

    global TargetSet_
    TargetSet_ = False

    global findMax_
    findMax_ = False

    # we also take the opporunity to set the random seed used by qbsolv. Qbsolv has a default random seed
    # so we mimic that behaviour here.
    if seed is None:
        seed = random.randint(0, 1L<<30)
    cdef int64_t c_seed = seed
    srand(c_seed)
    
    # ok, all of the globals are set, so let's get to actually solving the given problem. First we need the
    # list of variables used by Q
    variables = set().union(*Q)

    # we also set the inputs to qbsolv, using cython to make them proper C values
    cdef int n_variables = len(variables)
    cdef int n_repeats = repeats
    cdef int8_t **solution_list = <int8_t **>malloc2D(n_solutions + 1, n_variables, sizeof(int8_t))
    cdef double *energy_list = <double *>malloc((n_solutions + 1) * sizeof(double))
    cdef int *solution_counts = <int *>malloc((n_solutions + 1) * sizeof(int))
    cdef int *Qindex = <int *>malloc((n_solutions + 1) * sizeof(int))

    # create a matrix and set everything to 0
    cdef double **Qarr = <double **>malloc2D(n_variables, n_variables, sizeof(double))
    for row in range(n_variables):
        for col in range(n_variables):
            Qarr[row][col] = 0.

    cdef int u, v
    cdef double bias
    # put the values from Q into Qarr
    for (u, v), bias in iteritems(Q):
        # upper triangular
        if v < u:
            Qarr[v][u] = -bias
        else:
            Qarr[u][v] = -bias

    # Ok, solve using qbsolv! This puts the answer into output_sample
    solve(Qarr, n_variables, repeats, solution_list, energy_list, solution_counts, Qindex, n_solutions)

    # we have three things we are interested in, the samples, the energies and the number of times each
    # sample appeared
    samples = []
    energies = []
    counts = []

    # it is actually faster to use a while loop here and keep everything as a C object
    cdef int i = 0
    while i < n_solutions:
        soln_idx = Qindex[i]  # Qindex tracks where the order of the solutions

        # if no solutions were found then we can stop
        if solution_counts[soln_idx] == 0:
            break

        samples.append({v: int(solution_list[soln_idx][v]) for v in range(n_variables)})
        energies.append(float(energy_list[soln_idx]))
        counts.append(int(solution_counts[soln_idx]))

        i += 1

    # free the allocated variables
    free(solution_list)
    free(energy_list)
    free(solution_counts)
    free(Qindex);
    free(Qarr)

    return samples, energies, counts
