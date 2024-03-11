/**********************************************************************************
 *  afmisozerocore.c - library for SDW calculations at T=0
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
 ***********************************************************************************/

#include "afmisozerocore.h"
#include "nint.h"
#include "quartic.h" 



double occupation_wrapper( double x,
			   void * param )
{
  struct Parameters * p = param;
  p->T = x;
  return SUMMOR(heatkernel,p,SUMPREC);
}



double iso_SCE_wrapper( double x,
			void * param )
{
  struct Parameters * p = param;
  p->T = x;
  return SUMMOR(iso_SCE_heatkernel,p,SUMPREC);
}



double eigen_occ_wrapper( double x,
			  void * param )
{
  struct Parameters * p = param;
  p->T = x;
  return SUMMOR(eigen_occ_heatkernel,p,SUMPREC);
}



double
heatkernel ( double x,
	     double y,
	     double z,
	     void * param )
{
  /* T will hold the energy value omega */
  /* eta already has a placeholder */
  struct Parameters * p = param;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double A = ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk2q ) * ( p->T - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM) * ( .75*(p->T-nk)*(p->T-nkq) + (p->T-nk)*(p->T-nk3q) + .75*(p->T-nk2q)*(p->T-nk3q) )
    + .5625 * (p->J)*(p->J)*(p->J)*(p->J)*(p->magM)*(p->magM)*(p->magM)*(p->magM);
  double dAdw = ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk2q )
    + ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk3q )
    + ( p->T - nk ) * ( p->T - nk2q ) * ( p->T - nk3q )
    + ( p->T - nkq ) * ( p->T - nk2q ) * ( p->T - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM)*( 1.75*(2.*p->T-nk-nk3q) + .75*(2.*p->T-nkq-nk2q) );
  return exp( A*A/dAdw/dAdw/(-2.)/(p->eta)/(p->eta) ) / sqrt(2.*M_PI) / (p->eta) ;
}





double
mu_finder( void * par,
	   const double nc,
	   const double lbound,
	   const double ubound  )
{
  struct Parameters * p = par;
  double lower, upper;
  double target;
  double res;
  char * ecode = NULL;

#ifdef DEBUG
  fprintf(p->fout,":: %f\n",nc);
#endif
  lower = lbound;
  upper = ubound;
  do
    {
      target = (lower+upper)/2.;
      p->mu=target;
      res = quad(occupation_wrapper,p,-10.,target,0,REL_ERROR,&ecode);
      if ( ecode )
	{
#ifdef DEBUG
	  fprintf(p->fout,"Error code, %s.",ecode);
#endif
	  free(ecode);
	  ecode = NULL;
	}
#ifdef DEBUG
      fprintf(p->fout,":: (%1.6e,%1.6e) @ %1.8e => %1.8e\n",lower,upper,target,res);
      fflush(p->fout);
#endif
      if (NEAR(res,nc,PREC/10.))
	break;
      if (res < nc)
	{
	  lower=target;
	}
      else
	{
	  upper=target;
	}
    } while (NEAR(upper,lower,PREC)==0);
#ifdef DEBUG
  fprintf(p->fout,":::::: %1.8e\n",target);
#endif
  return target;
}




double
mu_stepper( void * par,
	    const double nc,
	    const double lbound,
	    const double ubound  )
{
  struct Parameters * p = par;
  double lower, upper;
  double inc, sum, testsum;
  double target;
  double res;
  char * ecode = NULL;
  
#ifdef DEBUG
  fprintf(p->fout,":: %f\n",nc);
#endif
  sum = 0.;
  lower = lbound;
  upper = ubound;
  inc = (upper-lower)/10.;
  target = lower + inc;
  while ( (NEAR(target,lower,PREC)==0) )
    {
#ifdef DEBUG
	  fprintf(p->fout,"[%f,%f] Sum=%1.6e, Testsum=%1.6e, Inc=%1.3e ",lower,target,sum,testsum,inc);
	  fflush(p->fout);
#endif
      testsum = sum;
      p->mu=target;
      testsum += quad(occupation_wrapper,p,lower,target,0,REL_ERROR,&ecode);
      if ( ecode )
	{
	  free(ecode);
	  ecode = NULL;
	}

      if ( NEAR(testsum,nc,PREC)==1 )
	{
	  break;
	}
      if (testsum < nc)
	{
	  lower  =  target;
	  target += inc;
	  sum = testsum;
	}
      else
	{
	  inc /= 2.;
	  target = lower + inc;
	}
      if ( target > upper )
	{
	  target = upper*2.;
	  break;
	}
    }
#ifdef DEBUG
  fprintf(p->fout,": %1.8e\n",target);
  fflush(p->fout);
#endif
  return target;
}






double
iso_SCE_zero( double x,
	      double y,
	      double z,
	      void * param )
{
  struct Parameters * p = param;
  int ret, i;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);
  double * r;
  double buff = 0.;
  double dAdw, denom;
  
  ret = quartic_solve(x,y,z,param,&r);
  if ( ret < 1 )
    {
      free(r);
      return 0.;
    }
  else
    {
      for ( i=0;i<ret;i++ )
	{
	  if ( p->mu > r[i] )
	    {
	      dAdw = ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q )
		+ ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk3q )
		+ ( r[i] - nk ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		+ ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		- jmsq*( 1.75*(2.*r[i]-nk-nk3q) + .75*(2.*r[i]-nkq-nk2q) );
	      denom = 1. -
		(  ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		   - jmsq * ( .75*(r[i]-nk)*(r[i]-nkq) + (r[i]-nk)*(r[i]-nk3q) + .75*(r[i]-nk2q)*(r[i]-nk3q) ) 
		   + .5625 * jmsq * jmsq ) * 
		( 2.*((r[i]-nk)*(r[i]-nk2q)-1.75*jmsq)+(2.*r[i]-nk-nk3q)*(2.*r[i]-nk2q-nkq)+2.*((r[i]-nk)*(r[i]-nk3q)-.75*jmsq)+(2.*r[i]-nkq-nk2q)*(2.*r[i]-nk3q-nk) ) /
		dAdw/dAdw;
	      buff += ( 1.5*(r[i]-nkq)*(r[i]-nk) + 2.*(r[i]-nk)*(r[i]-nk3q) + 1.5*(r[i]-nk2q)*(r[i]-nk3q) - 2.25*jmsq )
		/ (  (2.*r[i]-nk-nk3q) * ((r[i]-nkq)*(r[i]-nk2q)-1.75*jmsq) + (2.*r[i]-nkq-nk2q) * ((r[i]-nk)*(r[i]-nk3q)-.75*jmsq) ) /
		fabs(denom);
	    }
	}
    }
  free(r);
  return (p->J)*(p->magM)*buff; 
}



double
iso_SCE_heatkernel( double x,
		    double y,
		    double z,
		    void * param )
{
  struct Parameters * p = param;
  int ret, i;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);
  double B = 1.5*(p->T-nkq)*(p->T-nk) + 2.*(p->T-nk)*(p->T-nk3q) + 1.5*(p->T-nk2q)*(p->T-nk3q) - 2.25*jmsq ;
  double A = ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk2q ) * ( p->T - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM) * ( .75*(p->T-nk)*(p->T-nkq) + (p->T-nk)*(p->T-nk3q) + .75*(p->T-nk2q)*(p->T-nk3q) )
    + .5625 * (p->J)*(p->J)*(p->J)*(p->J)*(p->magM)*(p->magM)*(p->magM)*(p->magM);
  double dAdw = ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk2q )
    + ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk3q )
    + ( p->T - nk ) * ( p->T - nk2q ) * ( p->T - nk3q )
    + ( p->T - nkq ) * ( p->T - nk2q ) * ( p->T - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM)*( 1.75*(2.*p->T-nk-nk3q) + .75*(2.*p->T-nkq-nk2q) );
  return  B / dAdw * exp( A*A/dAdw/dAdw/(-2.)/(p->eta)/(p->eta) ) / sqrt(2.*M_PI) / (p->eta) ;
}



double
eigen_occ ( double x,
	    double y,
	    double z,
	    void * param )
{
  struct Parameters * p = param;
  int ret, i;
  double * r;
  double buff = 0.;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);
  double dAdw;
  double denom;

  ret = quartic_solve(x,y,z,param,&r);
  if ( ret < 1 )
    {
      free(r);
      return 0.;
    }
  else
    {
      for ( i=0;i<ret;i++ )
	{
	  if ( p->mu > r[i] )
	    {
	      dAdw = ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q )
		+ ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk3q )
		+ ( r[i] - nk ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		+ ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		- jmsq*( 1.75*(2.*r[i]-nk-nk3q) + .75*(2.*r[i]-nkq-nk2q) );
	      denom = 1. -
		(  ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		   - jmsq * ( .75*(r[i]-nk)*(r[i]-nkq) + (r[i]-nk)*(r[i]-nk3q) + .75*(r[i]-nk2q)*(r[i]-nk3q) ) 
		   + .5625 * jmsq * jmsq ) * 
		( 2.*((r[i]-nk)*(r[i]-nk2q)-1.75*jmsq)+(2.*r[i]-nk-nk3q)*(2.*r[i]-nk2q-nkq)+2.*((r[i]-nk)*(r[i]-nk3q)-.75*jmsq)+(2.*r[i]-nkq-nk2q)*(2.*r[i]-nk3q-nk) ) /
		dAdw/dAdw;
	      buff += r[i] / fabs(denom);
	    }
	}
    }
  free(r);
  return buff;
}



double
eigen_occ_heatkernel ( double x,
		       double y,
		       double z,
		       void * param )
{
  struct Parameters * p = param;
  int ret, i;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);

  double A = ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk2q ) * ( p->T - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM) * ( .75*(p->T-nk)*(p->T-nkq) + (p->T-nk)*(p->T-nk3q) + .75*(p->T-nk2q)*(p->T-nk3q) )
    + .5625 * (p->J)*(p->J)*(p->J)*(p->J)*(p->magM)*(p->magM)*(p->magM)*(p->magM);
  double dAdw = ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk2q )
    + ( p->T - nk ) * ( p->T - nkq ) * ( p->T - nk3q )
    + ( p->T - nk ) * ( p->T - nk2q ) * ( p->T - nk3q )
    + ( p->T - nkq ) * ( p->T - nk2q ) * ( p->T - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM)*( 1.75*(2.*p->T-nk-nk3q) + .75*(2.*p->T-nkq-nk2q) );
  return ( p->T ) * exp( A*A/dAdw/dAdw/(-2.)/(p->eta)/(p->eta) ) / sqrt(2.*M_PI) / (p->eta) ;

}



double rootkernel ( double x, double y, double z, void * param )
{
  struct Parameters * p = param;
  int ret, i,j;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);
  double * r;
  double buff = 0.;
  double dAdw ;
  double denom;
  
  ret = quartic_solve(x,y,z,param,&r);
  if ( ret < 1 )
    {
      free(r);
      return 0.;
    }
  else
    {
      for ( i=0;i<ret;i++ )
	{
	  if ( p->mu > r[i] )
	    {
	      denom = 1.;
	      for (j=0;j<ret;j++)
		{
		  if ( j != i )
		    {
		      denom *= r[i]-r[j] ;
		    }
		}
	      /* produces the same result as below */
	      dAdw = ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q )
		+ ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk3q )
		+ ( r[i] - nk ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		+ ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		- jmsq*( 1.75*(2.*r[i]-nk-nk3q) + .75*(2.*r[i]-nkq-nk2q) );
	      buff += dAdw / denom ;
	      
	      /* this attempt appears to be quite complicated, yet returns nothing spectacular:
		 If r[i] is a root of A, then 1 - A*(A'')/(A')^2 = 1 - 0 = 1.
		 And this is why all these calculations produce the same result as when I use buff+=1.
		 This implies that my iso_SCE function needs no modification, since its divisors
		 are also = 1.
		 But why, then, does THIS function not produce the correct nf ?
	      dAdw = ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q )
		+ ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk3q )
		+ ( r[i] - nk ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		+ ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		- jmsq*( 1.75*(2.*r[i]-nk-nk3q) + .75*(2.*r[i]-nkq-nk2q) );
	      denom = 1. -
		(  ( r[i] - nk ) * ( r[i] - nkq ) * ( r[i] - nk2q ) * ( r[i] - nk3q )
		   - jmsq * ( .75*(r[i]-nk)*(r[i]-nkq) + (r[i]-nk)*(r[i]-nk3q) + .75*(r[i]-nk2q)*(r[i]-nk3q) ) 
		   + .5625 * jmsq * jmsq ) * 
		( 2.*((r[i]-nk)*(r[i]-nk2q)-1.75*jmsq)+(2.*r[i]-nk-nk3q)*(2.*r[i]-nk2q-nkq)+2.*((r[i]-nk)*(r[i]-nk3q)-.75*jmsq)+(2.*r[i]-nkq-nk2q)*(2.*r[i]-nk3q-nk) ) /
		dAdw/dAdw;
		buff += fabs( 1./denom );*/

	      /* B = dAdw */
	    }
	}
    }
  free(r);
  return buff;
}
