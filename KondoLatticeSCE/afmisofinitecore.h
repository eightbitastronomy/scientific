/**********************************************************************************
 *  afmisofinitecore.h - Library for SDW phase calculations at T>0
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

#ifndef AFM_ISO_FINITE_H
#define AFM_ISO_FINITE_H


#include <stdlib.h>
#include "kdefs.h"
#include <math.h>
#include "nint.h"


enum { T,
       M,
       N,
       J,
       X,
       Y,
       Z  };



double occupation_c_boundary( double x,
			      double y,
			      double z,
			      void * param      );

double lindhard_boundary( double x,
			  double y,
			  double z,
			  void * param    );

double magRHS( void * param,
	       double (*f)(double,double,double,void *) );

double mufinder( struct Parameters * mp,
		 const double nc,
		 const double lbound,
		 const double ubound  );

int countq(int c);




#endif /*AFM_ISO_FINITE_H*/
