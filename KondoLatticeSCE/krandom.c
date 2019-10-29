/********************************************************************************************************************************
 *
 *  krandom library: library to provide pseudorandom functions; to hide access to mersenne twister from dSFMT library.
 *
 *  Author (pseudonomously): eightbitastronomy (eightbitastronomy@protonmail.com)
 *  KondoLatticeSCE is copyrighted by eightbitastronomy, 2018.
 *  dSFMT code is copyrighted by its owners.
 * 
 *  License information:
 *
 *  dSFMT is licensed separately from KondoLatticeSCE. See dSFMT/LICENSE.txt.
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
 ********************************************************************************************************************************
 *
 *  2018 aug 23: moved "shotgun" functions from aug23 version of kondo_bulk to here. Added "multiplier" argument to _slug.
 *
 *
 ********************************************************************************************************************************/





#include "krandom.h"
#include "dSFMT/dSFMT.h"
#include <time.h>
#include <math.h>


dsfmt_t mtwister;

void
start_krand()
{
  time_t seed;
  time( &seed );
  dsfmt_init_gen_rand(&mtwister,seed);
}


int
generate_q_partial( double * x,
		    double * y,
		    double * z )
{
  double * vector[3];
  vector[0] = x;
  vector[1] = y;
  vector[2] = z;
  int decision = dsfmt_genrand_uint32(&mtwister) % 3 ;
  int i;
  for ( i=0; i<=decision; i++ )
    {
      *(vector[2-i]) = M_PI * (  2.*dsfmt_genrand_close_open(&mtwister) - 1  );
    }
  return decision;
}


void
generate_q( double * x,
	    double * y,
	    double * z )
{
  *x = M_PI * (  2.*dsfmt_genrand_close_open(&mtwister) - 1  );
  *y = M_PI * (  2.*dsfmt_genrand_close_open(&mtwister) - 1  );
  *z = M_PI * (  2.*dsfmt_genrand_close_open(&mtwister) - 1  );
}


double
generate_unit()
{
  return dsfmt_genrand_close_open(&mtwister);
}



double
shotgun_slug( double a,
	      double b,
	      double multiplier )
{
  if ( a<b )
    {
      return ( b*dsfmt_genrand_close_open(&mtwister) + multiplier*b );
    }
  return ( a*dsfmt_genrand_close_open(&mtwister) + multiplier*a );
}



double
shotgun_shot( double a,
	      double b )
{
  if ( a<b )
    {
      return ( (b-a)*dsfmt_genrand_close_open(&mtwister) + a );
    }
  return ( (a-b)*dsfmt_genrand_close_open(&mtwister) + b );
}




void
stop_krand()
{
  /* nothing here as yet.
     dSFMT requires no destructor */
}
