/***********************************************************************
 *  This code is part of GLPK (GNU Linear Programming Kit).
 *
 *  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
 *  2009, 2010, 2011, 2013, 2017 Andrew Makhorin, Department for Applied
 *  Informatics, Moscow Aviation Institute, Moscow, Russia. All rights
 *  reserved. E-mail: <mao@gnu.org>.
 *
 *  GLPK is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  GLPK is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include "env.h"
#include "glpk.h"

#define SSXTRACE_DEAFULT_ITERATION_SPACE 8

glp_stmcp DEFAULT_STMCP = {
        .store_mem = GLP_STORE_TRACE_MEM_OFF,

        .basis_trace = GLP_BASIS_TRACE_OFF,
        .nonbasis_trace = GLP_NONBASIS_TRACE_OFF,
        .complexity_trace = GLP_COMPLEXITY_TRACE_ON,
        .objective_trace = GLP_OBJECTIVE_TRACE_ON,
        .pivot_rule = GLP_TRACE_PIVOT_DANTZIG,
        .bits_only = GLP_TRACE_BITS_ONLY_OFF,

        .info_file_basename = {'\0'},
        .objective_values_file_basename = {'\0'},
        .status_file_basename = {'\0'},
        .variable_values_file_basename = {'\0'},
};

glp_ssxtrace *glp_create_ssxtrace(const glp_stmcp *para) {
    glp_ssxtrace *info = glp_malloc(sizeof(glp_ssxtrace));

    if (para == NULL) {
        info->params = DEFAULT_STMCP;
    } else {
        info->params = *para;
    }

    info->no_basic = 0;
    info->no_nonbasic = 0;

    info->status = NULL;
    info->lb = NULL;
    info->ub = NULL;

    info->bases = NULL;
    info->basic_values = NULL;

    info->objective_values = NULL;

    info->no_iterations = 0;
    info->_allocated_iter = 0;

    info->updated = 0;

    info->info_fptr = NULL;
    if (strlen(para->info_file_basename) > 0) {
        info->info_fptr = fopen(para->info_file_basename, "w");
        xassert(info->info_fptr != NULL);
    }

    info->status_fptr = NULL;
    if (strlen(para->status_file_basename) > 0) {
        info->status_fptr = fopen(para->status_file_basename, "w");
        xassert(info->status_fptr != NULL);
    }

    info->objective_values_fptr = NULL;
    if (strlen(para->objective_values_file_basename) > 0) {
        info->objective_values_fptr = fopen(para->objective_values_file_basename, "w");
        xassert(info->objective_values_fptr != NULL);
    }

    info->variable_values_fptr = NULL;
    if (strlen(para->variable_values_file_basename) > 0) {
        info->variable_values_fptr = fopen(para->variable_values_file_basename, "w");
        xassert(info->variable_values_fptr != NULL);
    }

    return info;
}

void glp_ssxtrace_free(glp_ssxtrace *trace) {

    if (trace->status != NULL) { xfree(trace->status); }

    if (trace->lb != NULL) {
        size_t k = trace->no_basic + trace->no_nonbasic;
        for (size_t i = 1; i <= k; ++i) mpq_clear(trace->lb[i]);

        xfree(trace->lb);
    }

    if (trace->ub != NULL) {
        size_t k = trace->no_basic + trace->no_nonbasic;
        for (size_t i = 1; i <= k; ++i) mpq_clear(trace->ub[i]);

        xfree(trace->ub);
    }

    if (trace->bases != NULL) { xfree(trace->bases); }

    if (trace->basic_values != NULL) {

        // Number of allocated basic values
        size_t k = trace->no_iterations * trace->no_basic;
        for (size_t i = 0; i < k; ++i) mpq_clear(trace->basic_values[i]);

        xfree(trace->basic_values);
    }

    if (trace->objective_values != NULL) {
        for (size_t i = 0; i < trace->no_iterations; ++i)
            mpq_clear(trace->objective_values[i]);

        xfree(trace->objective_values);
    }

    if (trace->info_fptr != NULL)
        fclose(trace->info_fptr);

    if (trace->objective_values_fptr != NULL)
        fclose(trace->objective_values_fptr);

    if (trace->status_fptr != NULL)
        fclose(trace->status_fptr);

    if (trace->variable_values_fptr != NULL)
        fclose(trace->variable_values_fptr);
}

/* Ensure that there is enough space allocated for next iteration
 * of values. If there is not enough space, reallocate and double the current
 * size
 * */
void glp_ssxtrace_ensure_enough_space(glp_ssxtrace *trace) {
    /* There is enough space still */
    if (trace->no_iterations < trace->_allocated_iter) return;

    xassert(trace->no_basic > 0);
    xassert(trace->no_nonbasic > 0);

    if (trace->_allocated_iter == 0) {
        if (trace->params.nonbasis_trace == GLP_NONBASIS_TRACE_ON) {
            size_t k = trace->no_basic + trace->no_nonbasic + 1;
            trace->status =
                    xalloc(SSXTRACE_DEAFULT_ITERATION_SPACE * k, sizeof(int));
            xassert(trace->status != NULL);
        }

        if (trace->params.basis_trace == GLP_BASIS_TRACE_ON) {
            trace->bases =
                    xalloc(SSXTRACE_DEAFULT_ITERATION_SPACE * trace->no_basic,
                           sizeof(size_t));
            xassert(trace->bases != NULL);

            trace->basic_values =
                    xalloc(SSXTRACE_DEAFULT_ITERATION_SPACE * trace->no_basic,
                           sizeof(mpq_t));
            xassert(trace->basic_values != NULL);
        }

        if (trace->params.objective_trace == GLP_OBJECTIVE_TRACE_ON) {
            trace->objective_values =
                    xalloc(SSXTRACE_DEAFULT_ITERATION_SPACE, sizeof(mpq_t));
            xassert(trace->objective_values != NULL);
        }

        trace->_allocated_iter = SSXTRACE_DEAFULT_ITERATION_SPACE;

        return;
    }

    // Double the amount of space
    trace->_allocated_iter *= 2;

    if (trace->params.nonbasis_trace == GLP_NONBASIS_TRACE_ON) {
        size_t k = trace->no_basic + trace->no_nonbasic + 1;
        trace->status = xrealloc(trace->status, trace->_allocated_iter * k,
                                 sizeof(int));
        xassert(trace->status != NULL);
    }

    if (trace->params.basis_trace == GLP_BASIS_TRACE_ON) {
        trace->bases =
                xrealloc(trace->bases, trace->_allocated_iter * trace->no_basic,
                         sizeof(size_t));
        xassert(trace->bases != NULL);

        trace->basic_values = xrealloc(trace->basic_values,
                                       trace->_allocated_iter * trace->no_basic,
                                       sizeof(mpq_t));
        xassert(trace->basic_values != NULL);
    }

    if (trace->params.objective_trace == GLP_OBJECTIVE_TRACE_ON) {
        trace->objective_values = xrealloc(
                trace->objective_values, trace->_allocated_iter, sizeof(mpq_t));
        xassert(trace->objective_values != NULL);
    }
}

/* eof */