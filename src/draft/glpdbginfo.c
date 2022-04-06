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

#define DBGINFO_DEAFULT_ITERATION_SPACE 8

glp_dbginfo *glp_create_dbginfo(void) {
    glp_dbginfo *info = glp_malloc(sizeof(glp_dbginfo));

    info->no_basic = 0;
    info->no_nonbasic = 0;

    info->bases = NULL;
    info->basic_values = NULL;

    info->nonbases = NULL;
    info->nonbasic_values = NULL;

    info->objective_values = NULL;

    info->no_iterations = 0;
    info->_allocated_iter = 0;

    info->updated = 0;
    return info;
}

// TODO: Deadllocate mpq_t values
void glp_dbginfo_free(glp_dbginfo *info) {
    //xprintf("Freeing debug info\n");

    if (info->bases != NULL) {
        xfree(info->bases);
        //xprintf("Freed basis\n");
    }

    if (info->basic_values != NULL) {

        // Number of allocated basic values
        size_t k = info->no_iterations * info->no_basic;
        for (size_t i = 0; i < k; ++i)
            mpq_clear(info->basic_values[i]);

        xfree(info->basic_values);
        //xprintf("Freed basic values\n");
    }

    if (info->nonbases != NULL)
        xfree(info->nonbases);

    if (info->nonbasic_values != NULL) {
        // Number of allocated non-basic values
        //size_t k = info->no_iterations * info->no_nonbasic;

        //for (size_t i = 0; i < k; ++i)
        //    mpq_clear(info->basic_values[i]);

        xfree(info->nonbasic_values);
        //xprintf("Freed non-basic values\n");
    }

    if (info->objective_values != NULL) {
        for (size_t i = 0; i < info->no_iterations; ++i)
            mpq_clear(info->objective_values[i]);

        xfree(info->objective_values);

        //xprintf("Freed objective values\n");
    }

    //xprintf("Successfully freed\n");
}

/* Ensure that there is enough space allocated for next iteration
 * of values. If there is not enough space, reallocate and double the current
 * size
 * */
void glp_dbginfo_ensure_enough_space(glp_dbginfo *info) {
    /* There is enough space still */
    if (info->no_iterations < info->_allocated_iter)
        return;

    xassert(info->no_basic > 0);

    if (info->_allocated_iter == 0) {
        info->bases = xalloc(
                DBGINFO_DEAFULT_ITERATION_SPACE * info->no_basic,
                sizeof(size_t));
        xassert(info->bases != NULL);
        info->basic_values = xalloc(
                DBGINFO_DEAFULT_ITERATION_SPACE * info->no_basic,
                sizeof(mpq_t));
        xassert(info->basic_values != NULL);
        info->nonbases = xalloc(
                DBGINFO_DEAFULT_ITERATION_SPACE * info->no_nonbasic,
                sizeof(size_t));
        xassert(info->nonbases != NULL);
        info->nonbasic_values = xalloc(
                DBGINFO_DEAFULT_ITERATION_SPACE * info->no_nonbasic,
                sizeof(mpq_t));
        xassert(info->nonbasic_values != NULL);
        info->objective_values = xalloc(
                DBGINFO_DEAFULT_ITERATION_SPACE,
                sizeof(mpq_t));
        xassert(info->objective_values != NULL);

        info->_allocated_iter = DBGINFO_DEAFULT_ITERATION_SPACE;

        return;
    }

    // Double the amount of space
    info->_allocated_iter *= 2;

    info->bases = xrealloc(
            info->bases,
            info->_allocated_iter * info->no_basic,
            sizeof(size_t));
    xassert(info->bases != NULL);

    info->basic_values = xrealloc(
            info->basic_values,
            info->_allocated_iter * info->no_basic,
            sizeof(mpq_t));
    xassert(info->basic_values != NULL);

    info->nonbases = xrealloc(
            info->nonbases,
            info->_allocated_iter * info->no_nonbasic,
            sizeof(size_t));
    xassert(info->nonbases != NULL);

    info->nonbasic_values = xrealloc(
            info->nonbasic_values,
            info->_allocated_iter * info->no_nonbasic,
            sizeof(mpq_t));
    xassert(info->nonbasic_values != NULL);

    info->objective_values = xrealloc(
            info->objective_values,
            info->_allocated_iter,
            sizeof(mpq_t));
    xassert(info->objective_values != NULL);
}

/* eof */