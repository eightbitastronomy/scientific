/**********************************************************************************
 *  afmisozerocore.h - library for SDW calculations at T=0
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
#ifndef AFMISOZEROCORE_H
#define AFMISOZEROCORE_H





#include "kdefs.h"
#include <math.h>




double occupation_wrapper ( double x,
			    void * param );
double eigen_occ_wrapper ( double x,
			   void * param );
double iso_SCE_wrapper ( double x,
			 void * param );

double mu_finder( void * par,
		  const double nc,
		  const double lbound,
		  const double ubound );
double mu_stepper( void * par,
		   const double nc,
		   const double lbound,
		   const double ubound  );

double heatkernel ( double x,
		    double y,
		    double z,
		    void * param );
double rootkernel ( double x,
		    double y,
		    double z,
		    void * param );

double eigen_occ ( double x,
		   double y,
		   double z,
		   void * param );
double eigen_occ_heatkernel (
			     double x,
			     double y,
			     double z,
			     void * param );

double iso_SCE_zero( double x,
		     double y,
		     double z,
		     void * param );
double iso_SCE_heatkernel( double x,
			   double y,
			   double z,
			   void * param );








#endif /* AFMISOZEROCORE_H */
