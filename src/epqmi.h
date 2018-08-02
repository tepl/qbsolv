#ifndef EPQMI_H_INCLUDED
#define EPQMI_H_INCLUDED

/**********************************************************************
 *Proprietary Information D-Wave Systems Inc.                         *
 * Copyright (c) 2017 by D-Wave Systems Inc. All rights reserved.     *
 * Notice this code is licensed to authorized users only under the    *
 * applicable license agreement see eula.txt                          *
 * D-Wave Systems Inc., 3033 Beta Ave., Burnaby, BC, V5G 4M9, Canada. *
 **********************************************************************/

/*
 * This include file defines the interface to the epqmi library, which
 * is part of the qOp tools package.  Epqmi stands for "Embedded
 * Parameterized Quantum Machine Instruction".  The abstraction
 * supported by this library is a QUBO defined over named variables
 * and which contains optional parameters and assertions, and which
 * has been embedded into a D-Wave system with some specific geometry.
 *
 * Before using this library, create a QUBO using the format supported
 * by the dw(1) command and then embed the QUBO using the "dw embed"
 * subcommand.  The result is a binary workspace file with suffix
 * .eqpmi as well as a text workspace file with suffix .qubo, both of
 * which can be examined using dw subcommands.
 *
 * To use this library from a C program, set the dw workspace using
 * "dw set connection", "dw set solver", and "dw cd ..." before
 * starting the program.  Start the C program, which will use internal
 * dw environment variables to determine the current workspace.
 * Create an epqmi object using the DW_epqmi_read(...) function, which
 * returns a pointer to the opaque data type.
 *
 * All the rest of the library functions require a pointer to DW_epqmi
 * as the first argument.  These functions allow the C program to
 * query the underlying QUBO and find out its variables and
 * parameters.  The parameters can be bound to specified values to
 * create a QMI.  The QMI can be executed and characteristics of the
 * samples pulled from the QMI's distribution can then be examined.
 * When the epqmi object is no longer needed, free the memory
 * associated with it via DW_epqmi_free(...)
 */

/****************************************************************/

/*
 * The DW_epqmi data structure is opaque - meaning that all the
 * information associated with it is stored in fields which are
 * defined elsewhere.  The following definition includes only a header
 * field which is automatically populated with a magic value when the
 * epqmi is created.  Do not attempt to allocate a DW_epqmi.  The user
 * program should only contain pointers to DW_epqmi.  The only valid
 * way to allocate a DW_epqmi structure is with the DW_epqmi_read(...)
 * function.
 */

typedef struct DW_epqmi
{
  unsigned int magic;
} DW_epqmi;

/*
 * Use DW_epqmi_read(...) to create a DW_epqmi structure by reading
 * information from the current dw workspace.  If the epqmi_file
 * argument is a NULL pointer, assume the epqmi file is named
 * "default.epqmi".  Otherwise, assume that epqmi_file names an epqmi
 * file in the current workspace.  This returns NULL if there is a
 * problem and a non-NULL pointer if all is successful.
 */

DW_epqmi *DW_epqmi_read(char *epqmi_file);

/*
 * Use DW_epqmi_list_params(...) to list the parameters in the epqmi.
 * The param_names argument should be the address of a char **
 * variable.  The variable's value will be overwritten with the
 * address of an array of pointers to char.  The params argument
 * should be the address of an int which will be overwritten with the
 * number of parameters in the underlying QUBO.  At bind time, a value
 * must be specified for each parameter and the order of parameter
 * values must match the order of the parameter names returned by this
 * function.
 */

int DW_epqmi_list_params(DW_epqmi *epqmi, char ***param_names, int *params);

/*
 * Use DW_epqmi_list_vars(...) to list the variables in the epqmi.
 * The var_names argument should be the address of a char ** so that
 * the variable's value can be overwritten with the address of an
 * array of pointers to char.  The vars argument should be the
 * address of an int which will be overwritten with the number of
 * variables in the underlying QUBO.  After execution, the
 * DW_epqmi_sol_vars(...) function reports the value of each variable
 * in each sample.  The ordering of variables reported by that
 * function matches the ordering of variable names returned by this
 * function.
 */

int DW_epqmi_list_vars(DW_epqmi *epqmi, char ***var_names, int *vars);

/*
 * Use DW_epqmi_bind(...) to assign a value to each parameter in an
 * epqmi and create a QMI.  The parameter values are provided in an
 * array whose address is the second argument to the function.  The
 * order of parameter values in the array corresponds to the order of
 * parameter names returned by DW_epqmi_list_params(...).
 */

int DW_epqmi_bind(DW_epqmi *epqmi, float *param_values);

/*
 * Use DW_epqmi_exec(...) to execute the QMI created by
 * DW_epqmi_bind(...).  Specify the number of samples via the
 * num_reads argument.
 */

int DW_epqmi_exec(DW_epqmi *epqmi, int num_reads, int annealing_time);

/*
 * Use DW_epqmi_sols(...) after executing a QMI to determine the
 * number of unique samples in the distribution.  The solutions
 * argument is a pointer to a variable defined by the calling program
 * whose value will be overwritten with the number of distinct
 * samples.  If annealing_time is greater than 0, then it specifies
 * the annealing time in microseconds.
 */

int DW_epqmi_sols(DW_epqmi *epqmi, int *solutions);

/*
 * Use DW_epqmi_sol_vars(...) after executing a QMI to determine the
 * variable values in each variable in a specific solution.  Solnum is
 * the solution number.  Var is a pointer to an array allocated by the
 * calling program whose size in bytes equals the number of
 * variables.  This array will be overwritten with the 0 and 1 values
 * from the specified solution.  Additionally, the assertions
 * associated with the epqmi will be evaluated.  If all assertions are
 * true, the valid variable will be overwritten with 1, otherwise it
 * will be overwritten with 0.
 */

int DW_epqmi_sol_vars(DW_epqmi *epqmi, int solnum, char *var, int *valid);

/*
 * Use DW_epqmi_sol_occurs(...) after executing a QMI to determine how
 * many times a specific sample occurs.  Solnum names the solution and
 * the number of occurrences will be overwritten into the variable
 * whose address is passed in the occurrrences argument.
 */
int DW_epqmi_sol_occurs(DW_epqmi *epqmi, int solnum, int *occurrences);

/*
 * Use DW_epqmi_sol_obj(...) to get the objective value of a specific
 * sample.  The objective argument should be the address of a value
 * defined in the calling program which will be overwritten with the
 * objective of the specified solution.
 */
int DW_epqmi_sol_obj(DW_epqmi *epqmi, int solnum, float *objective);

/*
 * Use DW_epqmi_free(...) to free the storage associated with an
 * eqpqmi.
 */
int DW_epqmi_free(DW_epqmi *epqmi);

#endif // EPQMI_H_INCLUDED
