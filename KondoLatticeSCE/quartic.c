/**********************************************************************************
 *  quartic.c - library for quartic polynomial solution
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

#include "quartic.h"
#include "kdefs.h"
#include <math.h>



#include <gsl/gsl_poly.h>



int quartic_solve( double x,
		   double y,
		   double z,
		   void * param,
		   double ** roots )
{
  struct Parameters * p = param;
  double e0 = EN(x,y,z);
  double e1 = EN(x+p->qx,y+p->qy,z+p->qz);
  double e2 = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double e3 = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->magM)*(p->J)*(p->magM);
  double coef[4];
  double cubic[3];
  int retcubic,retquad;
  double quad[2];
  double rad1, rad2;
  int i;
  int nextroot;

  if ( roots == NULL )
    {
      return -1;
    }
  *roots = malloc(4*sizeof(double));
  if ( *roots == NULL )
    {
      return -2;
    }
  for ( i=0;i<4;i++ )
    {
      (*roots)[i] = 0.;
    }

  coef[0] = e0*e1*e2*e3 + jmsq*( .5625*jmsq - .75*e0*e1 - e0*e3 - .75*e2*e3 );
  coef[1] = -1. * ( e1*e2*e3 + e0*e2*e3 + e0*e1*e3 + e0*e1*e2 - jmsq*(1.75*e0+.75*e1+.75*e2+1.75*e3) );
  coef[2] = e0*e1 + e0*e2 + e0*e3 + e1*e2 + e1*e3 + e2*e3 - 2.5*jmsq ;
  coef[3] = -1. * ( e0 + e1 + e2 + e3 ) ;

  retcubic = gsl_poly_solve_cubic( -1.*coef[2],
				   coef[1]*coef[3]-4.*coef[0],
				   -1.*(coef[1]*coef[1]+coef[0]*coef[3]*coef[3]-4.*coef[0]*coef[2]),
				   &(cubic[0]),
				   &(cubic[1]),
				   &(cubic[2])  );

  if ( retcubic < 1 )
    {
      free(*roots);
      *roots = NULL;
      return 0;
    }

  nextroot = 0;
  
  for ( i=0;i<retcubic;i++ )
    {

      if ( nextroot > 3 )
	{
	  break;
	}
      
      rad1 = .25*coef[3]*coef[3] + cubic[i] - coef[2];
      rad2 = .25*cubic[i]*cubic[i] - coef[0];
      
      if ( (rad1<0.) || (rad2<0.) )
	{
	  /* Imaginary roots found for quartic */
	  continue;
	}
      
      retquad = gsl_poly_solve_quadratic( 1.,
					  .5*coef[3] - sqrt(rad1),
					  .5*cubic[i] - sqrt(rad2),
					  &(quad[0]),
					  &(quad[1])  );
      
      if ( retquad > 0 )
	{
	  if ( NEAR(quad[0]*quad[0]*quad[0]*quad[0]+coef[3]*quad[0]*quad[0]*quad[0]+coef[2]*quad[0]*quad[0]+coef[1]*quad[0]+coef[0],0.,PREC/10.) == 1 )
	    {
	      (*roots)[nextroot] = quad[0];
	      nextroot++;
	    }
	  if ( NEAR(quad[1]*quad[1]*quad[1]*quad[1]+coef[3]*quad[1]*quad[1]*quad[1]+coef[2]*quad[1]*quad[1]+coef[1]*quad[1]+coef[0],0.,PREC/10.) == 1 )
	    {
	      (*roots)[nextroot] = quad[1];
	      nextroot++;
	    }
	}
      
      
      retquad = gsl_poly_solve_quadratic( 1.,
					  .5*coef[3] + sqrt(rad1),
					  .5*cubic[i] + sqrt(rad2),
					  &(quad[0]),
					  &(quad[1])  );
      
      if ( retquad > 0 )
	{
	  if ( NEAR(quad[0]*quad[0]*quad[0]*quad[0]+coef[3]*quad[0]*quad[0]*quad[0]+coef[2]*quad[0]*quad[0]+coef[1]*quad[0]+coef[0],0.,PREC/10.) == 1 )
	    {
	      (*roots)[nextroot] = quad[0];
	      nextroot++;
	    }
	  if ( NEAR(quad[1]*quad[1]*quad[1]*quad[1]+coef[3]*quad[1]*quad[1]*quad[1]+coef[2]*quad[1]*quad[1]+coef[1]*quad[1]+coef[0],0.,PREC/10.) == 1 )
	    {
	      (*roots)[nextroot] = quad[1];
	      nextroot++;
	    }
	}
      
    }

  return nextroot;

}
