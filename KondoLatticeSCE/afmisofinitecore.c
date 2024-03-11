/**********************************************************************************
 *  afmisofinitecore.c - Library for SDW phase calculations at T>0
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


#include "afmisofinitecore.h"





double occupation_c_boundary( double x,
			      double y,
			      double z,
			      void * param   )
{
  struct Parameters * p = param;
  return (  FERMI( EN(x,y,z),p->mu,KB,p->T )
	    + FERMI( EN(x+p->qx,y+p->qy,z+p->qz),p->mu,KB,p->T )
	    + FERMI( EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz),p->mu,KB,p->T )
	    + FERMI( EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz),p->mu,KB,p->T )  ) ;
}



double lindhard_boundary( double x,
			  double y,
			  double z,
			  void * param    )
{
  struct Parameters * p = param;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  return 1.5*((FERMI(nk2q,p->mu,KB,p->T)-FERMI(nk3q,p->mu,KB,p->T))/(nk3q-nk2q))
    + 2.*((FERMI(nkq,p->mu,KB,p->T)-FERMI(nk2q,p->mu,KB,p->T))/(nk2q-nkq))
    + 1.5*((FERMI(nk,p->mu,KB,p->T)-FERMI(nkq,p->mu,KB,p->T))/(nkq-nk)) ;
}



double magRHS( void * param,
	       double (*f)(double,double,double,void *) )
{
  struct Parameters * p = param;
  return 0.9375*(p->J)*(p->J)*SUMMOR(f,p,SUMPREC)/KB/(p->T);
}



double mufinder( struct Parameters * mp,
		 const double nc,
		 const double lbound,
		 const double ubound  )
{
  double lower = lbound;
  double target;
  double res;
  double upper = ubound;

  do
    {
      target = (lower+upper)/2.;
      mp->mu=target;
      res = SUMMOR(occupation_c_boundary,mp,SUMPREC);
      if (res < nc)
	{
	  lower=target;
	}
      else
	{
	  upper=target;
	}
    } while (NEAR(upper,lower,PREC)==0);
  return target;
}



int countq(int c)
{
  return ((c*c*c)+(3*c*c)+(2*c))/6;
}
