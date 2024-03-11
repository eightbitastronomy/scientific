/*************************************************************************************************************************
 *  nint.c - Numerical computation & Encapsulating Library
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
 *  2018 aug 21: removed some comments
 *
 *****************************************************************************************************************************************/

#include "nint.h"




gsl_integration_cquad_workspace * w;
gsl_interp * wi;
gsl_interp_accel * acc;
double * xval, * yval;



/* n is # of integration intervals to give gsl a guide for the minimum necessary resources */
int start_nint( int n )
{
  w = gsl_integration_cquad_workspace_alloc((size_t)n);
  gsl_set_error_handler_off();
}


int start_ninterp( int m, double * x, double * y )
{
  wi = gsl_interp_alloc( gsl_interp_linear, (size_t)m );
  acc = gsl_interp_accel_alloc();
  xval = x;
  yval = y;
  return gsl_interp_init( wi, x, y, (size_t)m );
}


double mesh( double z )
{
  return gsl_interp_eval( wi, xval, yval, z, acc );
}




/* -------------------definition of quad:-----------------------------------------
 * int gsl_integration_cquad (const gsl_function * f,   the integrand
 *                            double a,                 lower limit
 *                            double b,                 upper limit
 *                            double epsabs,            abs error limit
 *                            double epsrel,            rel error limit
 *                            gsl_integration_cquad_workspace * workspace, 
 *                            double * result,          calculated result
 *                            double * abserr,          error in result
 *                            size_t * nevals   )       # evals used for result
 * ------------------------------------------------------------------------------- */

double quad(double (*f)(double,void *),void * p,double lb,double ub,double abs,double rel,char ** ecode)
{
  double ret;
  gsl_function F = { f, p };
  int retcode;
  double error;
  int nn;
  
  *ecode = NULL;
  retcode = gsl_integration_cquad(&F,lb,ub,abs,rel,w,&ret,NULL,NULL);
  if ( retcode )
    {
      *ecode = malloc(32*sizeof(char));
      /* for retcodes, see below */
      switch (retcode)
	{
	case GSL_EMAXITER:
	  strcpy(*ecode,"GSL: MAX ITERATIONS\0");
	  break;
	case GSL_EROUND:
	  strcpy(*ecode,"GSL: ROUNDING\0");
	  break;
	case GSL_ESING:
	  strcpy(*ecode,"GSL: SINGULARITY\0");
	  break;
	case GSL_EDIVERGE:
	  strcpy(*ecode,"GSL: DIVERGENT\0");
	  break;
	case GSL_EDOM:
	  strcpy(*ecode,"GSL: DOMAIN ERROR\0");
	  break;
	case GSL_ERANGE:
	  strcpy(*ecode,"GSL: RANGE ERROR\0");
	  break;
	case GSL_ENOMEM:
	  strcpy(*ecode,"GSL: NO MEMORY\0");
	  break;
	case GSL_EINVAL:
	  strcpy(*ecode,"GSL: INVALID ARG\0");
	  break;
	case GSL_EFAILED:
	  strcpy(*ecode,"GSL: GENERIC FAIL\0");
	  break;
	case GSL_ETOL:
	  strcpy(*ecode,"GSL: FAILED TO REACH TOLERANCE\0");
	  break;
	default:
	  sprintf(*ecode,"GSL: OTHER %d\0",retcode);
	  break;
	}
    }
  return ret;
}




void stop_nint( )
{
  if (w)
    gsl_integration_cquad_workspace_free(w);
}

void stop_ninterp( )
{
  if (wi)
    gsl_interp_free(wi);
  if (acc)
    gsl_interp_accel_free(acc);
  if (xval)
    {
      free(xval);
      xval=NULL;
    }
  if (yval)
    {
      free(yval);
      yval=NULL;
    }
}





/* GSL general error codes */
/*GSL_SUCCESS  = 0,
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,     iteration has not converged 
  GSL_EDOM     = 1,      input domain error, e.g sqrt(-1) 
  GSL_ERANGE   = 2,      output range error, e.g. exp(1e100) 
  GSL_EFAULT   = 3,      invalid pointer 
  GSL_EINVAL   = 4,      invalid argument supplied by user 
  GSL_EFAILED  = 5,      generic failure 
  GSL_EFACTOR  = 6,      factorization failed 
  GSL_ESANITY  = 7,      sanity check failed - shouldn't happen 
  GSL_ENOMEM   = 8,      malloc failed 
  GSL_EBADFUNC = 9,      problem with user-supplied function 
  GSL_ERUNAWAY = 10,     iterative process is out of control 
  GSL_EMAXITER = 11,     exceeded max number of iterations 
  GSL_EZERODIV = 12,     tried to divide by zero 
  GSL_EBADTOL  = 13,     user specified an invalid tolerance 
  GSL_ETOL     = 14,     failed to reach the specified tolerance 
  GSL_EUNDRFLW = 15,     underflow 
  GSL_EOVRFLW  = 16,     overflow  
  GSL_ELOSS    = 17,     loss of accuracy 
  GSL_EROUND   = 18,     failed because of roundoff error 
  GSL_EBADLEN  = 19,     matrix, vector lengths are not conformant 
  GSL_ENOTSQR  = 20,     matrix not square 
  GSL_ESING    = 21,     apparent singularity detected 
  GSL_EDIVERGE = 22,     integral or series is divergent 
  GSL_EUNSUP   = 23,     requested feature is not supported by the hardware 
  GSL_EUNIMPL  = 24,     requested feature not (yet) implemented 
  GSL_ECACHE   = 25,     cache limit exceeded 
  GSL_ETABLE   = 26,     table limit exceeded 
  GSL_ENOPROG  = 27,     iteration is not making progress towards solution 
  GSL_ENOPROGJ = 28,     jacobian evaluations are not improving the solution 
  GSL_ETOLF    = 29,     cannot reach the specified tolerance in F 
  GSL_ETOLX    = 30,     cannot reach the specified tolerance in X 
  GSL_ETOLG    = 31,     cannot reach the specified tolerance in gradient 
  GSL_EOF      = 32      end of file 
*/
  


double trapezoidal( double (*f)(double,double,double,void *),void * par,double targetinc )
{
  int n = ((int)(2.*M_PI/targetinc))+1;
  double inc = 2.*M_PI/(n-1);
  double v = inc*inc*inc;
  double sum=0.;
  double sub=0.;
  int i,j,k;
  /* reduced: */
  /* corners x 8 */
  sub += f( -M_PI, -M_PI, -M_PI, par );
  sub += f( -M_PI, -M_PI, -M_PI + ((double)(n-1)*inc), par );
  sub += f( -M_PI, -M_PI + ((double)(n-1)*inc), -M_PI, par );
  sub += f( -M_PI, -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI, -M_PI, par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI, -M_PI + ((double)(n-1)*inc), par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), -M_PI, par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), par );
  sum = .125*sub;
  /* edges x 12 */
  sub=0.;
  j=0; k=0;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0; k=0;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0; j=0;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=n-1; j=0;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  k=0;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  j=n-1; k=0;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  j=0; k=n-1;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=n-1; j=n-1;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=n-1; k=n-1;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  j=n-1; k=n-1;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  sum += .25*sub;
  /* faces x 6 */
  sub=0.;
  i=0;
  for (j=1;j<n-1;j++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  i=n-1;
  for (j=1;j<n-1;j++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  j=0;
  for (i=1;i<n-1;i++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  j=n-1;
  for (i=1;i<n-1;i++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  k=0;
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  k=n-1;
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  sum += .5*sub;
  /* bulk x 1 */
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  for (k=1;k<n-1;k++)
	    {
	      sum += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	    }
	}
	}
  return sum*v;
}


double k_trapezoidal( double (*f)(double,double,double,void *),void * par,double targetinc )
{
  int n = ((int)(2.*M_PI/targetinc))+1;
  double inc = 2.*M_PI/(n-1);
  double v = inc*inc*inc/8./M_PI/M_PI/M_PI;
  double sum=0.;
  double sub=0.;
  int i,j,k;
  /* reduced: */
  /* corners x 8 */
  sub += f( -M_PI, -M_PI, -M_PI, par );
  sub += f( -M_PI, -M_PI, -M_PI + ((double)(n-1)*inc), par );
  sub += f( -M_PI, -M_PI + ((double)(n-1)*inc), -M_PI, par );
  sub += f( -M_PI, -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI, -M_PI, par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI, -M_PI + ((double)(n-1)*inc), par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), -M_PI, par );
  sub += f( -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), -M_PI + ((double)(n-1)*inc), par );
  sum = .125*sub;
  /* edges x 12 */
  sub=0.;
  j=0; k=0;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0; k=0;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0; j=0;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=n-1; j=0;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  k=0;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  j=n-1; k=0;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  j=0; k=n-1;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=0;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=n-1; j=n-1;
  for (k=1;k<n-1;k++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  i=n-1; k=n-1;
  for (j=1;j<n-1;j++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  j=n-1; k=n-1;
  for (i=1;i<n-1;i++)
    {
      sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
    }
  sum += .25*sub;
  /* faces x 6 */
  sub=0.;
  i=0;
  for (j=1;j<n-1;j++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  i=n-1;
  for (j=1;j<n-1;j++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  j=0;
  for (i=1;i<n-1;i++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  j=n-1;
  for (i=1;i<n-1;i++)
    {
      for (k=1;k<n-1;k++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  k=0;
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  k=n-1;
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  sub += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	}
    }
  sum += .5*sub;
  /* bulk x 1 */
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  for (k=1;k<n-1;k++)
	    {
	      sum += f( -M_PI+((double)i*inc), -M_PI+((double)j*inc), -M_PI+((double)k*inc), par );
	    }
	}
	}
  return sum*v;
}



double lefthand_e( double (*f)(double,double,double,void *), double (*e)(double,double,double,void*),void * par,double inc,double * err )
{
  /* as a test-of-concept, this will be very specific. Range will be -Pi to Pi on all 3 variables.*/
  double v = inc*inc*inc;
  double kx=-M_PI;
  double ky=-M_PI;
  double kz=-M_PI;
  double sum=0.;
  double error;
  double retval;
  double upperbound=M_PI; /* if i use M_PI+inc, then the integral has the same value twice due to periodicity */
  do
    {
      ky=-M_PI;
      do
	{
	  kz=-M_PI;
	  do
	    {
	      sum+=f(kx,ky,kz,par);
	      retval=e(kx,ky,kz,par);
	      if ( retval > error )
		{
		  error = retval;
		}
	      kz+=inc;
	    } while (kz<upperbound);
	  ky+=inc;
	} while (ky<upperbound);
      kx+=inc;
    } while (kx<upperbound);
  *err = 3.*(2.*M_PI/inc)*v*error*error*error/12.;
  return sum*v;
}


double g_lefthand_c( double (*f)(double,double,double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double inc )
{
  double v = inc*inc*inc;
  double kx=lb;
  double ky;
  double kz;
  double sum=0.;
  do
    {
      ky=lb;
      do
	{
	  kz=lb;
	  do
	    {
	      sum+=f(kx,ky,kz,par);
	      kz+=inc;
	    } while (kz<ub);
	  ky+=inc;
	} while (ky<ub);
      kx+=inc;
    } while (kx<ub);
  return sum*v;
}



double trapezoidal_crude( double (*f)(double,double,double,void *),void * par,double inc,double * errorbuff )
{
  int n = ((int)(2.*M_PI/inc))+1;
  double v = inc*inc*inc;
  double sum=0.;
  double error;
  double max = 0.;
  double last,current,next;
  int i,j,k;
  /* brute-force calculation: */
  for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
	{
	  last = f( ((double)(i)*inc)-M_PI, ((double)(j)*inc)-M_PI, (-1.)*inc-M_PI, par );
	  current = f( ((double)(i)*inc)-M_PI, ((double)(j)*inc)-M_PI, -1.*M_PI, par );
	  for (k=1;k<n+1;k++)
	    {
	      next = f( ((double)(i)*inc)-M_PI, ((double)(j)*inc)-M_PI, ((double)k*inc)-M_PI, par );
	      sum += current;
	      error = (next - 2*current + last) / inc / inc;
	      if ( fabs(error) > fabs(max) )
		{
		  max = error;
		}
	      last = current;
	      current = next;
	    }
	}
    }
  *errorbuff = 3.*((double)n)*v*max*max*max/12. ;
  return v*sum;
}




/* this is actually "left-hand-q", but I have no reason to rewrite it" */

double trapezoidalq( __float128 (*f)(__float128,__float128,__float128,void *),void * par,__float128 inc )
{
  /* as a test-of-concept, this will be very specific. Range will be -Pi to Pi on all 3 variables.*/
  double v = inc*inc*inc;
  __float128 kx=-M_PI;
  __float128 ky=-M_PI;
  __float128 kz=-M_PI;
  __float128 sum=0.;
  __float128 upperbound=M_PI; /* if i use M_PI+inc, then the integral has the same value twice due to periodicity */
  do
    {
      ky=-M_PI;
      do
	{
	  kz=-M_PI;
	  do
	    {
	      sum+=f(kx,ky,kz,par);
	      kz+=inc;
	    } while (kz<upperbound);
	  ky+=inc;
	} while (ky<upperbound);
      kx+=inc;
    } while (kx<upperbound);
  return (double)((double)sum)*v;
}



/*********************************************
   Hasn't been checked. Use g_simpsons_c.
   -- v hasn't been updated,
   -- for loops haven't been checked...?
**********************************************/
double simpsons( double (*f)(double,double,double,void *),void * par,double targetinc )
{

  double sum=0.;
  double subsumA=0.;
  double subsumB=0.;
  double incbuff1,incbuff2,incbuff3;
  int i,j,k;
  int n = ((int)(2.*M_PI/targetinc))+1;
  double inc = 2.*M_PI/(n-1);
  double v = inc*inc*inc/216.;
  double half = inc/2.;
  

  /****** SURFACES ******/
  /* corners - are all vertex points */
  sum = f(-M_PI,-M_PI,-M_PI,par);
  sum += f(-M_PI,-M_PI,M_PI,par);
  sum += f(-M_PI,M_PI,-M_PI,par);
  sum += f(-M_PI,M_PI,M_PI,par);
  sum += f(M_PI,-M_PI,-M_PI,par);
  sum += f(M_PI,-M_PI,M_PI,par);
  sum += f(M_PI,M_PI,-M_PI,par);
  sum += f(M_PI,M_PI,M_PI,par);
  /* edges - (2 x vertex edges, 1 x half edges) */
  subsumB = f(-M_PI,-M_PI,half-M_PI,par);
  subsumB += f(-M_PI,M_PI,half-M_PI,par);
  subsumB += f(M_PI,-M_PI,half-M_PI,par);
  subsumB += f(M_PI,M_PI,half-M_PI,par);
  subsumB += f(-M_PI,half-M_PI,-M_PI,par);
  subsumB += f(-M_PI,half-M_PI,M_PI,par);
  subsumB += f(M_PI,half-M_PI,-M_PI,par);
  subsumB += f(M_PI,half-M_PI,M_PI,par);
  subsumB += f(half-M_PI,-M_PI,-M_PI,par);
  subsumB += f(half-M_PI,-M_PI,M_PI,par);
  subsumB += f(half-M_PI,M_PI,-M_PI,par);
  subsumB += f(half-M_PI,M_PI,M_PI,par);
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI;
      subsumA += f(-M_PI,-M_PI,incbuff1,par);
      subsumB += f(-M_PI,-M_PI,incbuff1+half,par);
      subsumA += f(-M_PI,M_PI,incbuff1,par);
      subsumB += f(-M_PI,M_PI,incbuff1+half,par);
      subsumA += f(M_PI,-M_PI,incbuff1,par);
      subsumB += f(M_PI,-M_PI,incbuff1+half,par);
      subsumA += f(M_PI,M_PI,incbuff1,par);
      subsumB += f(M_PI,M_PI,incbuff1+half,par);
      subsumA += f(-M_PI,incbuff1,-M_PI,par);
      subsumB += f(-M_PI,incbuff1+half,-M_PI,par);
      subsumA += f(-M_PI,incbuff1,M_PI,par);
      subsumB += f(-M_PI,incbuff1+half,M_PI,par);
      subsumA += f(M_PI,incbuff1,-M_PI,par);
      subsumB += f(M_PI,incbuff1+half,-M_PI,par);
      subsumA += f(M_PI,incbuff1,M_PI,par);
      subsumB += f(M_PI,incbuff1+half,M_PI,par);
      subsumA += f(incbuff1,-M_PI,-M_PI,par);
      subsumB += f(incbuff1+half,-M_PI,-M_PI,par);
      subsumA += f(incbuff1,-M_PI,M_PI,par);
      subsumB += f(incbuff1+half,-M_PI,M_PI,par);
      subsumA += f(incbuff1,M_PI,-M_PI,par);
      subsumB += f(incbuff1+half,M_PI,-M_PI,par);
      subsumA += f(incbuff1,M_PI,M_PI,par);
      subsumB += f(incbuff1+half,M_PI,M_PI,par);
    }
  sum += 2.*subsumA;
  sum += subsumB;
  /* faces - (4 x vertex faces, 2 x half-edges on face, 1 x half-faces on the face) */
  /* vertices on face */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  subsumA += f(-M_PI,incbuff1,incbuff2,par);
	  subsumA += f(M_PI,incbuff1,incbuff2,par);
	  subsumA += f(incbuff1,-M_PI,incbuff2,par);
	  subsumA += f(incbuff1,M_PI,incbuff2,par);
	  subsumA += f(incbuff1,incbuff2,-M_PI,par);
	  subsumA += f(incbuff1,incbuff2,M_PI,par);
	}
    }
  sum += 4.*subsumA;
  /* Half-edges on face */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    { 
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  subsumA += f(-M_PI,incbuff1,incbuff2,par);
	  subsumA += f(M_PI,incbuff1,incbuff2,par);
	  subsumA += f(-M_PI,incbuff2,incbuff1,par);
	  subsumA += f(M_PI,incbuff2,incbuff1,par);
	  subsumA += f(incbuff1,-M_PI,incbuff2,par);
	  subsumA += f(incbuff1,M_PI,incbuff2,par);
	  subsumA += f(incbuff2,-M_PI,incbuff1,par);
	  subsumA += f(incbuff2,M_PI,incbuff1,par);
	  subsumA += f(incbuff1,incbuff2,-M_PI,par);
	  subsumA += f(incbuff1,incbuff2,M_PI,par);
	  subsumA += f(incbuff2,incbuff1,-M_PI,par);
	  subsumA += f(incbuff2,incbuff1,M_PI,par);
	}
    }
  sum += 8.*subsumA; /* (x,y,1/2) on face gets multiplied by 4, and 2 for being on a face */
  /* Face-centers on face */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI+half;
	  subsumA += f(-M_PI,incbuff1,incbuff2,par);
	  subsumA += f(M_PI,incbuff1,incbuff2,par);
	  subsumA += f(-M_PI,incbuff2,incbuff1,par);
	  subsumA += f(M_PI,incbuff2,incbuff1,par);
	  subsumA += f(incbuff1,-M_PI,incbuff2,par);
	  subsumA += f(incbuff1,M_PI,incbuff2,par);
	  subsumA += f(incbuff2,-M_PI,incbuff1,par);
	  subsumA += f(incbuff2,M_PI,incbuff1,par);
	  subsumA += f(incbuff1,incbuff2,-M_PI,par);
	  subsumA += f(incbuff1,incbuff2,M_PI,par);
	  subsumA += f(incbuff2,incbuff1,-M_PI,par);
	  subsumA += f(incbuff2,incbuff1,M_PI,par);
	}
    }
  sum += 16.*subsumA; /* (x,1/2,1/2) on face gets multiplied by 16, and 1 for being on a face */  
  /****** BULK ******/
  /* [8 x vertex points
      1 x centers (1/2,1/2,1/2)
      2 x faces (x,1/2,1/2)
      4 x edges (x,y,1/2)       ]  */
  /* vertex points */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI,par);
	    }
	}
    }
  sum += 8.*subsumA;
  /* center points */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI+half;
	  for (k=0;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI+half,par);
	    }
	}
    }
  sum += 64.*subsumA; /* Bulk (1/2,1/2,1/2) get mult by 1 and by 64.*/
  /* bulk faces (x,1/2,1/2) */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI+half;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI,par);
	      subsumA += f(incbuff1,((double)k)*inc-M_PI,incbuff2,par);
	      subsumA += f(((double)k)*inc-M_PI,incbuff1,incbuff2,par);
	    }
	}
    }
  sum += 32.*subsumA; /* bulk (x,1/2,1/2) get mult by 2 and by 16 */
  /* bulk edges (x,y,1/2) */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI,par);
	      subsumA += f(incbuff1,((double)k)*inc-M_PI,incbuff2,par);
	      subsumA += f(((double)k)*inc-M_PI,incbuff1,incbuff2,par);
	    }
	}
    }
  sum += 16.*subsumA; /* bulk (x,y,1/2) get mult by 4 and by 4 */
  /****** finish ******/
  return sum*v;
}



double g_simpsons_c( double (*f)(double,double,double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double targetinc )
{
  double sum=0.;
  double subsumA=0.;
  double subsumB=0.;
  double incbuff1,incbuff2,incbuff3;
  int i,j,k;
  int n = ((int)((ub-lb)/targetinc))+1;
  double inc = (ub-lb)/(n-1);
  double v = inc*inc*inc/216.;
  double half = inc/2.;
  

  /****** SURFACES ******/
  /* corners - are all vertex points */
  sum = f(lb,lb,lb,par);
  sum += f(lb,lb,ub,par);
  sum += f(lb,ub,lb,par);
  sum += f(lb,ub,ub,par);
  sum += f(ub,lb,lb,par);
  sum += f(ub,lb,ub,par);
  sum += f(ub,ub,lb,par);
  sum += f(ub,ub,ub,par);
  /* edges - (2 x vertex edges, 1 x half edges) */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1=((double)i)*inc+lb+half;
      subsumA += f(lb,lb,incbuff1,par);
      subsumA += f(lb,ub,incbuff1,par);
      subsumA += f(ub,lb,incbuff1,par);
      subsumA += f(ub,ub,incbuff1,par);
      subsumA += f(lb,incbuff1,lb,par);
      subsumA += f(lb,incbuff1,ub,par);
      subsumA += f(ub,incbuff1,lb,par);
      subsumA += f(ub,incbuff1,ub,par);
      subsumA += f(incbuff1,lb,lb,par);
      subsumA += f(incbuff1,lb,ub,par);
      subsumA += f(incbuff1,ub,lb,par);
      subsumA += f(incbuff1,ub,ub,par);
    }
  sum += 4.*subsumA;
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1=((double)i)*inc+lb;
      subsumA += f(lb,lb,incbuff1,par);
      subsumA += f(lb,ub,incbuff1,par);
      subsumA += f(ub,lb,incbuff1,par);
      subsumA += f(ub,ub,incbuff1,par);
      subsumA += f(lb,incbuff1,lb,par);
      subsumA += f(lb,incbuff1,ub,par);
      subsumA += f(ub,incbuff1,lb,par);
      subsumA += f(ub,incbuff1,ub,par);
      subsumA += f(incbuff1,lb,lb,par);
      subsumA += f(incbuff1,lb,ub,par);
      subsumA += f(incbuff1,ub,lb,par);
      subsumA += f(incbuff1,ub,ub,par);
    }
  sum += 2.*subsumA;
  /* faces - (4 x vertex faces, 2 x half-edges on face, 1 x half-faces on the face) */
  /* vertices on face */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc+lb;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb;
	  subsumA += f(lb,incbuff1,incbuff2,par);
	  subsumA += f(ub,incbuff1,incbuff2,par);
	  subsumA += f(incbuff1,lb,incbuff2,par);
	  subsumA += f(incbuff1,ub,incbuff2,par);
	  subsumA += f(incbuff1,incbuff2,lb,par);
	  subsumA += f(incbuff1,incbuff2,ub,par);
	}
    }
  sum += 4.*subsumA;
  /* Half-edges on face */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    { 
      incbuff1 = ((double)i)*inc+lb+half;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb;
	  subsumA += f(lb,incbuff1,incbuff2,par);
	  subsumA += f(ub,incbuff1,incbuff2,par);
	  subsumA += f(lb,incbuff2,incbuff1,par);
	  subsumA += f(ub,incbuff2,incbuff1,par);
	  subsumA += f(incbuff1,lb,incbuff2,par);
	  subsumA += f(incbuff1,ub,incbuff2,par);
	  subsumA += f(incbuff2,lb,incbuff1,par);
	  subsumA += f(incbuff2,ub,incbuff1,par);
	  subsumA += f(incbuff1,incbuff2,lb,par);
	  subsumA += f(incbuff1,incbuff2,ub,par);
	  subsumA += f(incbuff2,incbuff1,lb,par);
	  subsumA += f(incbuff2,incbuff1,ub,par);
	}
    }
  sum += 8.*subsumA; /* (x,y,1/2) on face gets multiplied by 4, and 2 for being on a face */
  /* Face-centers on face */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc+lb+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb+half;
	  subsumA += f(lb,incbuff1,incbuff2,par);
	  subsumA += f(ub,incbuff1,incbuff2,par);
	  subsumA += f(incbuff1,lb,incbuff2,par);
	  subsumA += f(incbuff1,ub,incbuff2,par);
	  subsumA += f(incbuff1,incbuff2,lb,par);
	  subsumA += f(incbuff1,incbuff2,ub,par);
	}
    }
  sum += 16.*subsumA; /* (x,1/2,1/2) on face gets multiplied by 16, and 1 for being on a face */  
  /****** BULK ******/
  /* [8 x vertex points
      1 x centers (1/2,1/2,1/2)
      2 x faces (x,1/2,1/2)
      4 x edges (x,y,1/2)       ]  */
  /* vertex points */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc+lb;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc+lb,par);
	    }
	}
    }
  sum += 8.*subsumA;
  /* center points */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc+lb+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb+half;
	  for (k=0;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc+lb+half,par);
	    }
	}
    }
  sum += 64.*subsumA; /* Bulk (1/2,1/2,1/2) get mult by 1 and by 64.*/
  /* bulk center-faces (x,1/2,1/2) */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc+lb+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb+half;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc+lb,par);
	      subsumA += f(incbuff1,((double)k)*inc+lb,incbuff2,par);
	      subsumA += f(((double)k)*inc+lb,incbuff1,incbuff2,par);
	    }
	}
    }
  sum += 32.*subsumA; /* bulk (x,1/2,1/2) get mult by 2 and by 16 */
  /* bulk edges (x,y,1/2) */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc+lb;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc+lb;
	  for (k=0;k<n-1;k++)
	    {
	      incbuff3 = ((double)k)*inc+lb+half;
	      subsumA += f(incbuff1,incbuff2,incbuff3,par);
	      subsumA += f(incbuff1,incbuff3,incbuff2,par);
	      subsumA += f(incbuff3,incbuff1,incbuff2,par);
	    }
	}
    }
  sum += 16.*subsumA; /* bulk (x,y,1/2) get mult by 4 and by 4 */
  /****** finish ******/
  return sum*v;
}


double k_simpsons_c( double (*f)(double,double,double,void *),
		     void * par,
		     double targetinc )
{
  double sum=0.;
  double subsumA=0.;
  double subsumB=0.;
  double incbuff1,incbuff2,incbuff3;
  int i,j,k;
  int n = ((int)((2.*M_PI)/targetinc))+1;
  double inc = (2.*M_PI)/(n-1);
  double v = inc*inc*inc/1728/M_PI/M_PI/M_PI; /* inc^3 / (216*8pi^3) = inc^3 / 1728pi^3 */
  double half = inc/2.;
  
  /****** SURFACES ******/
  /* corners - are all vertex points */
  sum = f(-M_PI,-M_PI,-M_PI,par);
  sum += f(-M_PI,-M_PI,M_PI,par);
  sum += f(-M_PI,M_PI,-M_PI,par);
  sum += f(-M_PI,M_PI,M_PI,par);
  sum += f(M_PI,-M_PI,-M_PI,par);
  sum += f(M_PI,-M_PI,M_PI,par);
  sum += f(M_PI,M_PI,-M_PI,par);
  sum += f(M_PI,M_PI,M_PI,par);
  /* edges - (2 x vertex edges, 1 x half edges) */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1=((double)i)*inc-M_PI+half;
      subsumA += f(-M_PI,-M_PI,incbuff1,par);
      subsumA += f(-M_PI,M_PI,incbuff1,par);
      subsumA += f(M_PI,-M_PI,incbuff1,par);
      subsumA += f(M_PI,M_PI,incbuff1,par);
      subsumA += f(-M_PI,incbuff1,-M_PI,par);
      subsumA += f(-M_PI,incbuff1,M_PI,par);
      subsumA += f(M_PI,incbuff1,-M_PI,par);
      subsumA += f(M_PI,incbuff1,M_PI,par);
      subsumA += f(incbuff1,-M_PI,-M_PI,par);
      subsumA += f(incbuff1,-M_PI,M_PI,par);
      subsumA += f(incbuff1,M_PI,-M_PI,par);
      subsumA += f(incbuff1,M_PI,M_PI,par);
    }
  sum += 4.*subsumA;
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1=((double)i)*inc-M_PI;
      subsumA += f(-M_PI,-M_PI,incbuff1,par);
      subsumA += f(-M_PI,M_PI,incbuff1,par);
      subsumA += f(M_PI,-M_PI,incbuff1,par);
      subsumA += f(M_PI,M_PI,incbuff1,par);
      subsumA += f(-M_PI,incbuff1,-M_PI,par);
      subsumA += f(-M_PI,incbuff1,M_PI,par);
      subsumA += f(M_PI,incbuff1,-M_PI,par);
      subsumA += f(M_PI,incbuff1,M_PI,par);
      subsumA += f(incbuff1,-M_PI,-M_PI,par);
      subsumA += f(incbuff1,-M_PI,M_PI,par);
      subsumA += f(incbuff1,M_PI,-M_PI,par);
      subsumA += f(incbuff1,M_PI,M_PI,par);
    }
  sum += 2.*subsumA;
  /* faces - (4 x vertex faces, 2 x half-edges on face, 1 x half-faces on the face) */
  /* vertices on face */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  subsumA += f(-M_PI,incbuff1,incbuff2,par);
	  subsumA += f(M_PI,incbuff1,incbuff2,par);
	  subsumA += f(incbuff1,-M_PI,incbuff2,par);
	  subsumA += f(incbuff1,M_PI,incbuff2,par);
	  subsumA += f(incbuff1,incbuff2,-M_PI,par);
	  subsumA += f(incbuff1,incbuff2,M_PI,par);
	}
    }
  sum += 4.*subsumA;
  /* Half-edges on face */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    { 
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  subsumA += f(-M_PI,incbuff1,incbuff2,par);
	  subsumA += f(M_PI,incbuff1,incbuff2,par);
	  subsumA += f(-M_PI,incbuff2,incbuff1,par);
	  subsumA += f(M_PI,incbuff2,incbuff1,par);
	  subsumA += f(incbuff1,-M_PI,incbuff2,par);
	  subsumA += f(incbuff1,M_PI,incbuff2,par);
	  subsumA += f(incbuff2,-M_PI,incbuff1,par);
	  subsumA += f(incbuff2,M_PI,incbuff1,par);
	  subsumA += f(incbuff1,incbuff2,-M_PI,par);
	  subsumA += f(incbuff1,incbuff2,M_PI,par);
	  subsumA += f(incbuff2,incbuff1,-M_PI,par);
	  subsumA += f(incbuff2,incbuff1,M_PI,par);
	}
    }
  sum += 8.*subsumA; /* (x,y,1/2) on face gets multiplied by 4, and 2 for being on a face */
  /* Face-centers on face */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI+half;
	  subsumA += f(-M_PI,incbuff1,incbuff2,par);
	  subsumA += f(M_PI,incbuff1,incbuff2,par);
	  subsumA += f(incbuff1,-M_PI,incbuff2,par);
	  subsumA += f(incbuff1,M_PI,incbuff2,par);
	  subsumA += f(incbuff1,incbuff2,-M_PI,par);
	  subsumA += f(incbuff1,incbuff2,M_PI,par);
	}
    }
  sum += 16.*subsumA; /* (x,1/2,1/2) on face gets multiplied by 16, and 1 for being on a face */  
  /****** BULK ******/
  /* [8 x vertex points
      1 x centers (1/2,1/2,1/2)
      2 x faces (x,1/2,1/2)
      4 x edges (x,y,1/2)       ]  */
  /* vertex points */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI,par);
	    }
	}
    }
  sum += 8.*subsumA;
  /* center points */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI+half;
	  for (k=0;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI+half,par);
	    }
	}
    }
  sum += 64.*subsumA; /* Bulk (1/2,1/2,1/2) get mult by 1 and by 64.*/
  /* bulk center-faces (x,1/2,1/2) */
  subsumA = 0.;
  for (i=0;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI+half;
      for (j=0;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI+half;
	  for (k=1;k<n-1;k++)
	    {
	      subsumA += f(incbuff1,incbuff2,((double)k)*inc-M_PI,par);
	      subsumA += f(incbuff1,((double)k)*inc-M_PI,incbuff2,par);
	      subsumA += f(((double)k)*inc-M_PI,incbuff1,incbuff2,par);
	    }
	}
    }
  sum += 32.*subsumA; /* bulk (x,1/2,1/2) get mult by 2 and by 16 */
  /* bulk edges (x,y,1/2) */
  subsumA = 0.;
  for (i=1;i<n-1;i++)
    {
      incbuff1 = ((double)i)*inc-M_PI;
      for (j=1;j<n-1;j++)
	{
	  incbuff2 = ((double)j)*inc-M_PI;
	  for (k=0;k<n-1;k++)
	    {
	      incbuff3 = ((double)k)*inc-M_PI+half;
	      subsumA += f(incbuff1,incbuff2,incbuff3,par);
	      subsumA += f(incbuff1,incbuff3,incbuff2,par);
	      subsumA += f(incbuff3,incbuff1,incbuff2,par);
	    }
	}
    }
  sum += 16.*subsumA; /* bulk (x,y,1/2) get mult by 4 and by 4 */
  /****** finish ******/
  return sum*v;
}
  


double g_simpsons_l( double (*f)(double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double targetinc )
{
  double sum,subsum;
  int i;
  int n = ((int)((ub-lb)/targetinc))+1;
  double inc = (ub-lb)/(n-1);
  double v = inc/6.;
  double half = inc/2.;
  sum = f(lb,par) + f(ub,par);
  /* interval ends / vertices */
  subsum=0.;
  for (i=1;i<n-1;i++)
    {
      subsum += f(lb+((double)i)*inc,par);
    }
  sum += 2.*subsum;
  /* midpoints */
  subsum = 0.;
  for (i=0;i<n-1;i++)
    {
      subsum += f(lb+((double)i)*inc+half,par) ;
    }
  sum += 4.*subsum;
  return sum*v;
}


double g_simpsons_s( double (*f)(double,double,void *),
		     void * par,
		     double lb,
		     double ub,
		     double targetinc  )
{
  double sum,subsum;
  int i,j;
  int n = ((int)((ub-lb)/targetinc))+1;
  double inc = (ub-lb)/(n-1);
  double v = inc*inc/36./(ub-lb);
  double half = inc/2.;
  /* corners */
  sum = f(lb,lb,par) + f(lb,ub,par) + f(ub,lb,par) + f(ub,ub,par);
  /* edge, interval vertices */
  subsum=0.;
  for (i=1;i<n-1;i++)
    {
      subsum += f(lb+((double)i)*inc,lb,par);
      subsum += f(lb+((double)i)*inc,ub,par);
      subsum += f(lb,lb+((double)i)*inc,par);
      subsum += f(ub,lb+((double)i)*inc,par);
    }
  sum += 2.*subsum;
  /* edge, halfs */
  subsum=0.;
  for (i=0;i<n-1;i++)
    {
      subsum += f(lb+((double)i)*inc+half,lb,par);
      subsum += f(lb+((double)i)*inc+half,ub,par);
      subsum += f(lb,lb+((double)i)*inc+half,par);
      subsum += f(ub,lb+((double)i)*inc+half,par);
    }
  sum += 4.*subsum; /* 1x4xval */
  /* face, vertices */
  subsum = 0.;
  for (i=1;i<n-1;i++)
    {
      for (j=1;j<n-1;j++)
	{
	  subsum += f(lb+((double)i)*inc,lb+((double)j)*inc,par) ;
	}
    }
  sum += 4.*subsum;
  /* face, half-edges */
  subsum=0.;
  for (i=1;i<n-1;i++)
    {
      for (j=0;j<n-1;j++)
	{
	  subsum += f(lb+((double)i)*inc,lb+((double)j)*inc+half,par) ;
	  subsum += f(lb+((double)j)*inc+half,lb+((double)i)*inc,par) ;
	}
    }
  sum += 8.*subsum; /* 2x4xval */
  /* face, mid-face */
  subsum = 0.;
  for (i=0;i<n-1;i++)
    {
      for (j=0;j<n-1;j++)
	{
	  subsum += f(lb+((double)i)*inc+half,lb+((double)j)*inc+half,par) ;
	}
    }
  sum += 16.*subsum; /* 1x16xval */
  return sum*v;
}
