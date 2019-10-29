/**********************************************************************************
 *  bulkcore.h - Library for phase bulk/nonboundary calculations
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
 *******************************************************************************************************************************
 *
 *  2018 aug 20: attempts to correct Kondo energy have failed; implementing free_energy_kondo2 as a remedy.
 *  2018 aug 15: earlier changes to afm free energy fctns were not reflected here. Corrected.
 *
 *********************************************************************************************************************/

#ifndef BULK_H
#define BULK_H





#include <math.h>
#include "kdefs.h"
#include <stdlib.h>





double occupation_kondo_c( double x,
			   double y,
			   double z,
			   void * param );
double occupation_afm_c_wrapper( double x,
				 void * param );



double occupation_kondo_f( double x,
			   double y,
			   double z,
			   void * param );
double occupation_afm_f( double g,
			 double mu,
			 double j,
			 double s,
			 double t );



double sce_kondo( double x,
		  double y,
		  double z,
		  void * param );
double sce_afm_m( double g,
		  double mu,
		  double j,
		  double s,
		  double t );
double sce_afm_s_wrapper( double x,
			  void * param );
double sce_afm_s( double x,
		  double y,
		  double z,
		  void * param );




double heatkernel ( double x,
		    double y,
		    double z,
		    void * param );
double heatkernel_reduced ( double x,
			    double y,
			    double z,
			    void * param );



double free_energy_kondo ( double x,
			   double y,
			   double z,
			   void * param );
double free_energy_kondo2 ( double x,
			    double y,
			    double z,
			    void * param );
double free_energy_afm_band_wrapper ( double x,
				      void * param );
double free_energy_afm( void * param );




double mu_finder( double (*mf)(void *,double),
		  void * par,
		  double targetnc );
double mu_wrapper_afm( void * par,
		       double trn );
double mu_wrapper_kondo( void * par,
			 double trn );
double glam_finder( double (*gf)(void *, double),
		    void * par,
		    double nimp );
double glam_wrapper_afm( void * par,
			 double trg );
double glam_wrapper_kondo( void * par,
			   double trg );




int filter_minimum( double * record,
		    int recordsz,
		    int fields,
		    int key             );



  



#endif /*BULK_H*/
