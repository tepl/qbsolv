Qbsolv (no tabu search)
=======================

This version of qbsolv has no Tabu search and is supposed to be used
in DW environment only. Thus, one can see how good are solutions obtained
solely using hardware, i.e. without any help from classic solvers.

In addition, this version has couple more options to give a user
more control over the algorithm. Original qbsolv has these options hard-coded.

Added options:
 * [-p] paramChain
 * [-R] numReads
 * [-a] annealingTime

Removed options:
 * [-a] alrorithm
 * [-T] target mode
 * [-t] timeout
 * [-S] subproblemSize


Updated content of original README.txt
<pre>


This directory contains the following files and directories:
- README.txt:       this file
- LICENSE:          A copy of the Apache License Version 2.0
- src:              Source directory of the qbsolv program with Makefile
            that creates a binary read and solve "QUBO" files
- doc:          OpenOffice, PDF, and HTML versions of the man page,
            and instructions for how to change them
- example:      Directory of example(s) application(s) using qbsolv
            as a solver
- tests:        Directory of scripts and qubo(s) to test qbsolv
- contrib.txt:      Instructions for potential contributors


    qbsolv -i infile [-o outfile] [-m] [-n repeats] [-w] 
        [-h] [-v verbosityLevel] [-V] [-q] [-r seed]
        [-p paramChain] [-R numReads] [-a annealingTime]

DESCRIPTION 
    qbsolv executes a quadratic unconstrained binary optimization 
    (QUBO) problem represented in a file, providing bit-vector 
    result(s) that minimizes (or optionally, maximizes) the value of 
    the objective function represented by the QUBO.  The problem is 
    represented in the QUBO(5) file format and notably is not limited 
    to the size or connectivity pattern of the D-Wave system on which 
    it will be executed. 
        The options are as follows: 
    -i infile 
        The name of the file in which the input QUBO resides.  This 
        is a required option. 
    -o outfile 
        This optional argument denotes the name of the file to 
        which the output will be written.  The default is the 
        standard output. 
    -m 
        This optional argument denotes to find the maximum instead 
        of the minimum. 
    -n repeats 
        This optional argument denotes, once a new optimal value is 
        found, to repeat the main loop of the algorithm this number
        of times before stopping. The default value is 50.
    -w 
        If present, this optional argument will print the QUBO 
        matrix and result in .csv format. 
    -h 
        If present, this optional argument will print the help or 
        usage message for qbsolv and exit without execution. 
    -v verbosityLevel 
        This optional argument denotes the verbosity of output. A 
        verbosityLevel of 0 (the default) will output the number of 
        bits in the solution, the solution, and the energy of the 
        solution.  A verbosityLevel of 1 will output the same 
        information for multiple solutions, if found. A 
        verbosityLevel of 2 will also output more detailed 
        information at each step of the algorithm. This increases   
        the output up to a value of 4.
    -V 
        If present, this optional argument will emit the version 
        number of the qbsolv program and exit without execution. 
    -q 
        If present, this optional argument triggers printing the 
        format of the QUBO file.
    -r seed 
        Used to reset the seed for the random number generation.
    -p paramChain 
        This optional argument denotes strength of qubit chains 
        in a complete graph. The default value is 15.
    -R numReads 
        This optional argument denotes number of reads.
        The default value is 10.
    -a annealingTime 
        This optional argument denotes annealing time in microseconds.
        The default value is 20.


------------------------
qbsolv "qubo" input file format

   A .qubo file contains data which describes an unconstrained
quadratic binary optimization problem.  It is an ASCII file comprised
of four types of lines:

1) Comments - defined by a "c" in column 1.  They may appear
    anywhere in the file, and are otherwise ignored.

2) One program line, which starts with p in the first column.
    The program line must be the first non-comment line in the file.
    The program line has six required fields (separated by space(s)),
    as in this example:

  p   qubo  topology   maxNodes   nNodes   nCouplers

    where:
  p         the problem line sentinel
  qubo      identifies the file type
  topology  a string which identifies the topology of the problem
            and the specific problem type.  For an unconstrained problem,
            target will be "0" or "unconstrained."   Possible, for future
            implementations, valid strings might include "chimera128"
            or "chimera512" (among others).
  maxNodes   number of nodes in the topology.
  nNodes     number of nodes in the problem (nNodes <= maxNodes).
            Each node has a unique number and must take a value in the
            the range {0 - (maxNodes-1)}.  A duplicate node number is an
            error.  The node numbers need not be in order, and they need
            not be contiguous.
  nCouplers  number of couplers in the problem.  Each coupler is a
            unique connection between two different nodes.  The maximum
            number of couplers is (nNodes)^2.  A duplicate coupler is
            an error.

3) nNodes clauses.  Each clause is made up of three numbers.  The
            numbers are separated by one or more blanks.  The first two
            numbers must be integers and are the number for this node
            (repeated).  The node number must be in {0 , (maxNodes-1)}.
            The third value is the weight associated with the node, may be
            an integer or float, and can take on any positive or negative
            value, or zero.

4) nCouplers clauses.  Each clause is made up of three numbers.  The
            numbers are separated by one or more blanks.  The first two
            numbers must be different integers and are the node numbers
            for this coupler.  The two values (i and j) must have (i < j).
            Each number must be one of the nNodes valid node numbers (and
            thus in {0, (maxNodes-1)}).  The third value is the strength
            associated with the coupler, may be an integer or float, and can
            take on any positive or negative value, but not zero.  Every node
            must connect with at least one other node (thus must have at least
            one coupler connected to it).

Here is a simple QUBO file example for an unconstrained QUBO with 4
nodes and 6 couplers.  This example is provided to illustrate the
elements of a QUBO benchmark file, not to represent a real problem.

        | <--- column 1
        c
        c  This is a sample .qubo file
        c  with 4 nodes and 6 couplers
        c
        p  qubo  0  4  4  6 
        c ------------------
        0  0   3.4
        1  1   4.5
        2  2   2.1
        3  3   -2.4
        c ------------------
        0  1   2.2
        0  2   3.4
        1  2   4.5
        0  3   -2
        0  2   3.4
        1  2   4.5
        0  3   -2
        1  3   4.5678
        2  3   -3.22

</pre>

