/**********************************************************************************
 *  bare.c - library for bare interaction computations.
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
#include "bare.h"


extern double mesh(double);
extern double quad(double (*f)(double,void *), void *, double, double, double, double, char **);




/**************************************** Bare DoS & energy *******************************************/




double
dirac(double x)
{
  if ( x<=-6.0 )
    return 0.;
  if ( x>=6.0 )
    return 0.;
  return mesh(x);
}


double
diracw(double x,
       void * param)
{
  if ( x<=-6.0 )
    return 0.;
  if ( x>=6.0 )
    return 0.;
  return mesh(x);
}




/****************************** Kinetic Energies without mu-solutions ***********************************/




double
kinetic_bare_zero( double x,
		   void * param )
{
  return dirac(x)*x;
}


double
kinetic_bare_finite( double x,
		     void * param )
{
  struct Parameters * p = param;
  return FERMI(x,p->mu,KB,p->T)*dirac(x)*x;
}


double
kinetic_bare_zero_k(double x,
		    double y,
		    double z,
		    void * param )
{
  struct Parameters * p = param;
  double nrg = EN(x,y,z);
  if ( p->mu > nrg )
    return nrg;
  return 0.;
}


double
kinetic_bare_finite_k(double x,
		      double y,
		      double z,
		      void * param )
{
  struct Parameters * p = param;
  double nrg = EN(x,y,z);
  return FERMI(nrg,p->mu,KB,p->T)*nrg;
}


double
kinetic_bare_zero_k_q(double x,
		      double y,
		      double z,
		      void * param )
{
  struct Parameters * p = param;
  double nrg = EN(x+p->qx,y+p->qy,z+p->qz);
  if ( p->mu > nrg )
    return nrg;
  return 0.;
}


double
kinetic_bare_finite_k_q(double x,
			double y,
			double z,
			void * param )
{
  struct Parameters * p = param;
  double nrg = EN(x+p->qx,y+p->qy,z+p->qz);
  return FERMI(nrg,p->mu,KB,p->T)*nrg;
}




/*************************************** Paramagnetic Energy for Helmholtz Calculation **********************************************/




double
paramagnetic_energy( void   * par,
		     double lb,
		     double ub,
		     double * nc ,
		     double * mu    )
{
  struct Parameters * p = par;
  double retdbl;
  double target;
  double lower = lb;
  double upper = ub;
  char * ecode = NULL;
  int retint;

  /* find mu given no interaction and given nc as passed in */

  do
    {
      target = (lower+upper)/2.;
      retdbl = quad(occupation_bare_zero,NULL,lb,target,0,REL_ERROR,&ecode);
      if (retdbl < p->nc)
	{
	  lower=target;
	}
      else
	{
	  upper=target;
	}
    } while (NEAR(upper,lower,PREC)==0);

  /* with mu, calculate the PM energy */

  if ( nc )
    *nc = retdbl;
  if ( mu )
    *mu = target;

  retdbl = quad(kinetic_bare_zero,NULL,lb,target,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
    }

  return retdbl;  

}





/************************************* occupation number functions ********************************************/




double
occupation_bare_zero( double x,
		      void * par )
{
  return dirac(x);
}


double
occupation_bare_zero_k(double x,
		       double y,
		       double z,
		       void * par)
{
  struct Parameters * p = par;
  double nk = EN(x,y,z);
  if ( p->mu > nk )
    return 1.;
  else
    return 0.;
}


double
occupation_bare_finite( double x,
			void * par )
{
  struct Parameters * p = par;
  return FERMI(x,p->mu,KB,p->T)*dirac(x);
}


double
occupation_bare_finite_k( double x,
			  double y,
			  double z,
			  void * par )
{
  struct Parameters * p = par;
  return FERMI(EN(x,y,z),p->mu,KB,p->T);
}

