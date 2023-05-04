# glpk-simplex-trace

This repository is a fork of the GLPK (GNU Linear Programming Kit) Version 4.64
library Refer to [GLPK webpage](http://www.gnu.org/software/glpk/glpk.html).
It has been extended with the functionality to trace the basic solutions of the
exact simplex algorithm provided with the original distribution.

## Installation

The library uses CMake, so you can include it to a project using the standard
way. We assume the library source can be found in the
``libs/glpk-simplex-trace`` directory.

1. Include the library's source directory in your main ``CMakeLists.txt`` file:

```cmake
add_subdirectory(libs/glpk-simplex-trace/)
```

This tells CMake to process the ``CMakeLists.txt`` file located in the
``libs/glpk-simplex-trace/`` directory.

2. Link the library to your target:

Assuming you have already defined a target (e.g., an executable) called
``my_target``, link the ``glpk-simplex-trace`` library to it using the
``target_link_libraries()`` command:

```cmake
target_link_libraries(my_target PRIVATE glpk-simplex-trace)
```

The ``PRIVATE`` keyword indicates that the ``glpk-simplex-trace`` library is
only needed for building and linking the ``my_target`` target, but not for
other targets that may depend on it.

## Usage

The library builds on the standard GLPK API, please refer to the GLPK manual in
``doc/glpk.pdf`` for more details. The main method running the exact solver
with trace is the ``glp_exact_trace`` method. This method takes as argument the
``glp_prob`` object containing the problem, ``glp_smcp`` object of standard
control parameters for GLPK (refer to API) and ``glp_ssxtrace`` tracer object containing.
The tracer is then initialized using the ``glp_create_ssxtrace`` method, by supplying it with
``glp_stmcp`` object containing tracer control parameters.

The possible tracer control parameters are:

- ``store_mem``: Controls whether to store the solution in memory (Phase II **only**).

  Possible values are:
    * ``GLP_STORE_TRACE_MEM_ON`` -- save solution in memory (default)
    * ``GLP_STORE_TRACE_MEM_OFF``  -- do not save solution in memory (default)

- ``objective_trace``: Controls whether to store the objective values.

  Possible values are:
    * ``GLP_OBJECTIVE_TRACE_ON``: Turns tracing of objective values ON. (default)
    * ``GLP_OBJECTIVE_TRACE_OFF``: Turns tracing of objective values OFF.

- ``basis_trace``: Controls whether to store the values of basic variables.

  Possible values are:
    * ``GLP_BASIS_TRACE_ON``: Turns tracing of basic variable values ON. (default)
    * ``GLP_BASIS_TRACE_OFF``: Turns tracing of basic variable values OFF.

- ``status_trace``: Controls whether to store the status of variables.

  Possible values are:
    * ``GLP_STATUS_TRACE_ON``: Turns tracing of basic variable values ON. (default)
    * ``GLP_STATUS_TRACE_OFF``: Turns tracing of basic variable values OFF.

- ``pivot_rule``: Pivoting rule to use

  Possible values are:
    * ``GLP_TRACE_PIVOT_BEST``: To use the best (largest increase) rule.
    * ``GLP_TRACE_PIVOT_BLAND``: To use the Bland's rule.
    * ``GLP_TRACE_PIVOT_DANTZIG``: To use the Dantzig's (textbook) rule. (default)
    * ``GLP_TRACE_PIVOT_RANDOM``: To use the random rule (chooses entering and leaving variable at random).

- ``fractionality_bits_trace``: Whether to only story fractionality bits (FOR INTERNAL PURPOSES)

  Possible values are:
    * ``GLP_TRACE_BITS_ONLY_ON``: Turn bit tracing ON.
    * ``GLP_TRACE_BITS_ONLY_OFF``: Turn bit tracing OFF. (default)

- ``scale``: Scale the problem by multiplying RHS by the LCM of denominators (FOR INTERNAL PURPOSES).

  Possible values are:
    * ``GLP_TRACE_SCALE_ON``: Turn scaling ON.
    * ``GLP_TRACE_SCALE_OFF``: Turn scaling OFF. (default)

- ``info_file_basename``, ``objective_values_file_basename``, ``status_file_basename``, ``variable_values_file_basename``:
  Name of file where to store values. If left empty, the values will not be stored. All values empty by default.

NOTE: We do not recommend using the ``store_memory`` option, as the results can
get very large very quickly even for relatively small problems. It only
supports the phase II and storing the solution to a file is better in almost
all cases.

Usage example:

```c
void solve(char* filename) {
    glp_prob *P;
    P = glp_create_prob();

    // Use default Simplex Method Control Parameters (SMCP)
    glp_smcp params = NULL;

    // Read the problem from a file in an LP format
    glp_read_lp(P, NULL, filename);

    // Construct initial basis
    glp_adv_basis(P, 0);

    // Use Simplex Trace Method Control Parameters to default
    glp_stmcp trace_params = GLP_DEFAULT_STMCP;

    // Create a trace object for storing information
    glp_ssxtrace* trace = glp_create_ssxtrace(&trace_params);

    // Solve using trace
    glp_exact_trace(P, &params, trace);

    // Cleanup
    glp_ssxtrace_free(trace);
    glp_delete_prob(P);
}
```

For a complete example refer to ``dkubek/glpsol_tracer``

## Output Format

In this section we describe the format of the ouput files. The output format of a rational number
is either a integer or two integers separated by a slash ``/``

Status of a variable is one of 5 possible values

* ``1`` - basic variable
* ``2`` - non-basic variable on lower bound
* ``3`` - non-basic variable on upper bound
* ``4`` - non-basic free (unbounded) variable
* ``5`` - non-basic fixed variable

Variable can be either _logical_ or _structural_ (for more information on structural and logical variables refore to
GLPK
Manual)

- **Info file ouput**

  The info file output is split into two sections. The first section is right at the beginning of the file, with each
  line containing a ``key : value`` pair. The possible values listed are
    * ``rows``: The number of rows of the problem matrix. Also represents the number of logical variables.
    * ``cols``: The number of columns of the problem matrix. Also represents the number of structural variables.
    * ``nonzeros``: The number of nonzero values in the problem matrix.
    * ``iterations``: The number of iterations needed to solve the problem.
    * ``status``: The exit status of the solution.
    * ``objective``: The final objective value as a rational.
    * ``scale`` (optional): Integer by which the problem has been scaled (if the ``--scale`` option has been used).
      Dividing the objective value by this number yields the objective value of the original problem.

  The second part is delimited by
    ```
    --- BEGIN VARIABLES (TYPE-NAME-STATUS-VALUE) ---
    ...
    --- END VARIABLES ---
    ```
  Each line in this section contains three values separated by space in this order -- type, name, status, value
    * ``type``: denotes whether the variable is logical ``l`` or structural ``s``
    * ``name``: denotes the name of the variable
    * ``status``: denotes the final status of the variable
    * ``value``: denotes the final value of the variable

- **Variable file output format**

  The i-th line of the variable file contains a list of rational numbers separated by space representing the values of
  the variables in the i-th iteration. The order of variables is as defined in the info file.

- **Objective file output**
  The i-th line contains a single rational number representing the objective value in i-th iteration.

- **Status file output**
  The i-th line of the variable file contains a list of integers separated by space representing the status of
  the variables in the i-th iteration. The order of variables is as defined in the info file.