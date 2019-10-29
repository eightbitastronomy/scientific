/**********************************************************************************
 *  kfinitecore.h - library for Kondo calculations at T>0
 *
 *
 *  Author (pseudonomously): eightbitastronomy (eightbitastronomy@protonmail.com)
 *  Copyrighted by eightbitastronomy, 2018.
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

#ifndef KFINITECORE_H

#define KFINITECORE_H




#include <math.h>
#include <stdlib.h>
#include "kdefs.h"
#include "bare.h"



/* hybridization self-consistency equation, returns 1/J */
double sce( double x,
	    void * params );
double sce_k( double x,
	      double y,
	      double z,
	      void * params );
double sceV0( double x,
	      void * params );
double sceV0_k( double x,
		double y,
		double z,
		void * params );
double sceV0_k_eta( double x,
		double y,
		double z,
		void * params );

/* Fermi-Dirac function */
double fermi( const double xx,
	      const double mm,
	      const double TT   );

double sumd( double * a, int n );


double free_energy( double x,
		    void * params );


#endif /* KFINITECORE_H */
