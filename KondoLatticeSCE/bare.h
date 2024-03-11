/**********************************************************************************
 *  bare.h - library for bare interaction computations.
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
#ifndef BARE_H

#define BARE_H





#include "kdefs.h"
#include <math.h>
#include <stdlib.h>


/* implementation of dirac delta function */
double dirac(double x);

/* wrapper for dirac() for use with direct numerical integration */
double diracw(double x,
	      void * param);

/* conduction electron kinetic energy: epsilon_k c-dagger c */
double kinetic_bare_zero( double x,
			  void * param );
double kinetic_bare_finite( double x,
			    void * param );
double kinetic_bare_zero_k(double x,
			   double y,
			   double z,
			   void * param );
double kinetic_bare_finite_k(double x,
			     double y,
			     double z,
			     void * param );
double kinetic_bare_zero_k_q(double x,
			     double y,
			     double z,
			     void * param );
double kinetic_bare_finite_k_q(double x,
			       double y,
			       double z,
			       void * param );

/* paramagnetic helmholtz energy */
double paramagnetic_energy( void   * par,
			    double lb,
			    double ub,
			    double * nc ,
			    double * mu    );

/* conduction electron occupation numbers. No zero-k appears because it isn't practicable and isn't necessary */
double occupation_bare_zero( double x,
			     void * par );
double occupation_bare_zero_k(double x,
			      double y,
			      double z,
			      void * par);
double occupation_bare_finite( double x,
			       void * par );
double occupation_bare_finite_k( double x,
				 double y,
				 double z,
				 void * par );



#endif /*BARE_H*/
