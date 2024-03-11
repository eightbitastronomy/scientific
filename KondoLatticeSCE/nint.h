/*************************************************************************************************************************
 *  nint.h - Numerical computation & Encapsulating Library
 *
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
 *****************************************************************************************************************************
 *
 *  Encapsulating Library for 'calc' main.c
 *  The purpose of this library is simply to hide the GNU GSL mechanics from main.c
 *  For interpolation of the d.o.s. mesh, 
 *    1. call start_ninterp
 *    2. call mesh whenever an interpolated value is required
 *    3. call stop_ninterp when the interpolation is no longer needed.
 *  For the quadrature-based numerical integration,
 *    1. call start_nint
 *    2. call quad whenever an integration is required
 *    3. call stop_nint when no more integration is needed.
 *  The underlying GSL calls for mesh and quad require different function types and parameters, 
 *    so be sure to look at nint.c and the GSL documentation.
 *  Although GSL documentation says standard CBLAS libraries should work, only linking with "gslcblas" has succeeded
 *  1 oct 2016
 *  updates, 7 july 2017, added quadmath version, trapezoidalq
 ************************************************************************************************************************/


#ifndef NINTEGRATION


#define NINTEGRATION

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_errno.h"

int start_nint( int n );
int start_ninterp( int m,
		   double * x,
		   double * y );

double mesh( double z );
double quad(double (*f)(double,void *),
	    void * p,
	    double lb,
	    double ub,
	    double abs,
	    double rel,
	    char ** ecode);

void stop_nint( );
void stop_ninterp( );

/* 3D trapezoidal rule integration: 
      minimal grid calculations, no error calculations. */
double trapezoidal( double (*f)(double,double,double,void *),
		    void * par,
		    double inc );
double k_trapezoidal( double (*f)(double,double,double,void *),
		      void * par,
		      double inc );
/* 3D trapezoidal rule integration: 
      brute-force grid calculations, error bound computed from sampling in z direction only */
double trapezoidal_crude( double (*f)(double,double,double,void *),
			  void * par,
			  double inc,
			  double * errorbuff );
/* formerly trapezoidal: gives same value as trapezoidal, so i modified it for error bounds */
/*   The error function found in sdw_b_nqq.c, passed via e(), gave error bounds O(10^23) or so. Unacceptable.
     However, because I used an analytical expression for the error (2nd derivative of the Lindhard susc), possibly my expression is incorrect.
     Or, possibly, the max error lies on a delta-spike, which throws off the whole error calculation. */
double lefthand_e( double (*f)(double,double,double,void *),
		   double (*e)(double,double,double,void *), /* error calculator */
		   void * par,
		   double inc,
		   double * err );                           /* must be non-null and valid */
/* crude integration with lefthand rule over a cube on the interval [lb,ub]x[lb,ub]x[lb,ub]. */
double g_lefthand_c( double (*f)(double,double,double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double inc );

/* should be 'lefthandq' */
/*   It seems that on the tested CPU, a xeon, intel has only implemented quad precision in software, not hardware.
     The calculation slowed down tremendously (I aborted a short run after ~ 39 hours), leading me to conclude that
     either the use of quad-precision precludes 4-core parallelism or the accuracy gained in this method is not worth
     the resource cost. */
double trapezoidalq( __float128 (*f)(__float128,__float128,__float128,void *),
		     void * par,
		     __float128 inc );

/* simpson's 1/3 rule for 3D, tailored for a [-Pi,Pi] cube domain */
/* targetinc is used to calculate a practicable increment size which allows intervals to actually span [-pi,pi] */
double simpsons( double (*f)(double,double,double,void *),
		 void * par,
		 double targetinc  );
/* simpson's 1/3 rule for 3D: generalized bounds on a cube of [lb,ub]x[lb,ub]x[lb,ub] */
/* targetinc is used to calculate a practicable increment size which allows intervals to actually span [-pi,pi] */
double g_simpsons_c( double (*f)(double,double,double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double targetinc  );
double k_simpsons_c( double (*f)(double,double,double,void *),
		     void * par,
		     double targetinc );
double g_simpsons_l( double (*f)(double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double targetinc  );
double g_simpsons_s( double (*f)(double,double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double targetinc  );



#endif
