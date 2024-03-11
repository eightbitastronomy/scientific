/**********************************************************************************
 *  sm.h -- Parallel-processing library for multicore cpus.
 *          Uses forked processes and shared memory.
 *          Dispatch functions with input and output arrays.
 *          Set the number of logical cores here; it will not be autodetected.
 *
 *  Author: Miguel Abele
 *  Copyrighted by Miguel Abele, 2018.
 * 
 *  License information: 
 *
 *  This file is a part of KondoLatticeSCE.
 *
 *  KondoLatticeSCE is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 3
 *  of the License, or (at your option) any later version.
 *
 *  KondoLatticeSCE is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 ***********************************************************************************/

#ifndef SM_H
#define SM_H


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include <sys/unistd.h>
#include <sys/wait.h>
#include "list.h"


#define MAXTASK     3
#define SLEEPTIME   10000

struct d_Function {
  void * (*f)( void *, void *, void * ); /* params, input data, return data */
  void * param;
};


/* returns > 0 for child processes.
   Is tailored to kpms_rhs.          */
int dispatch( double (*f)(void *,double,double,double), void * params, double * data, int n, double ** result );


/* extension/generalization, this function uses a callback function so that children can exit immediately upon completion 
   and uses pointers for arbitrary data handling. However, it is still not general purpose, and is limited to zkpms.        */
void dispatchv( void * (*f)(void *,void *), void * params, void * data, int n, int sz, void * result, void (*callback)(int,...));


/* extension of dispatchv to general input and result data sizes, as well as handling of SIGQUIT (but not SIGINT) */
void dispatchv2( struct d_Function dfunc,          /* struct for function to be called on data as well as parameters */
		 int n,                            /* number of data points to be processed */
		 void * dataptr,                   /* the data points */
		 int dsz,                          /* size of each data point */
		 void * resultptr,                 /* storage of results */
		 int rsz,                          /* size of each result; allows for multiple quantities to be calculated from one set of data-point + parameters */
		 void (*callback)(int,...)    );   /* the exit function, expedites clean exits and proper handling of children-processes */


#endif /* SM_H */
