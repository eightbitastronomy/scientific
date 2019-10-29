#ifndef KZEROCORE_H

#define KZEROCORE_H





#include "kdefs.h"
#include "bare.h"
#include <math.h>
#include <stdlib.h>



/* impurity electron density of states */
double rho_m(double x,
	     void * param);
double rho_m_k(double x,     
	       double y,
	       double z,
	       void * params );


/* conduction electron density of states */
double rho_c(double x,
	     void * param);
double rho_c_k(double x,     
	       double y,
	       double z,
	       void * params );


/* hybridized particle density of states - MUST BE VERIFIED */
double rho_hyb( double x,      
		void * param);


/* Free energy hybridized terms: */
double eigen_occup( double x,
		    void * param );
double eigen_occup_k( double x,
		      double y,
		      double z,
		      void * param );

/* sceV() with explicit implementation of Heavisides instead of macros */
double sceV(double x,   /*J32 MUST BE CORRECTED 28NOV2017 */
	     void * p);
double sce_k(double x,  
	     double y,
	     double z,
	     void * param);

/* to be used as function handed to sm.c's dispatch routine */
void * sceVcaller(void * param,
		  void *dat,
		  void *retdat);

/* mu-solving: first stage, find flat region of occ vs mu. call msolve accordingly. */
int msplit(struct Parameters * par,
	   double lb,
	   double ub,
	   double * lmu,
	   double * umu);
int msplit_k(struct Parameters * par,
	     double lb,
	     double ub,
	     double * lmu,
	     double * umu);

/* mu-solving: second stage, bisecting search for mu. */
double msolve( struct Parameters * par,
	       const double lb,
	       const double ub,
	       const double offset );  /* shift test values: value += offset*eta, eta is defined as a fraction of PREC */
double msolve_k( struct Parameters * par,
		 const double lb,
		 const double ub,
		 const double offset );  /* shift test values: value += offset*eta, eta is defined as a fraction of PREC */

/* mu-solving: mu==gamma bisecting-search solver */
int mgsolve( struct Parameters * par,
	     double lb,
	     double ub );



#endif /* KZEROCORE_H */
