#include "kzerocore.h"


extern double quad(double (*f)(double,void *), void *, double, double, double, double, char **);
extern double SUMMOR(double (*f)(double,double,double,void *),void *,double);




/****************************************** Dressed DoS *********************************************/




double
rho_m(double x,     
      void * param)
{
  struct Parameters * p = param;
#ifdef USECASEJ32
  return (14.0625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g)/(x-p->g))*dirac(x - (14.0625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g))) ;
#else /*S=1/2*/
  return (.5625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g)/(x-p->g))*dirac(x - (.5625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g))) ;
#endif /*USECASEJ32*/
}


double
rho_c(double x,       
      void * param)
{
  struct Parameters * p = param;
#ifdef USECASEJ32
  return dirac(x - (14.0625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g))) ;
#else /*S=1/2*/
  return dirac(x - (.5625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g))) ;
#endif /*USECASEJ32*/
}


double rho_c_k(double x,       /* I haven't made the effort to consolidate the S=1/2 & S=3/2 if-then constructs into one...*/
	       double y,
	       double z,
	       void * params )
{
  struct Parameters * p = params;
  double nrg = EN(x,y,z);
  double buff = 0.;
#ifdef USECASEJ32
  double rad = sqrt( .25*(p->g-nrg)*(p->g-nrg)+14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  if ( p->mu > .5*(p->g+nrg)+rad )
    {
      buff += fabs(.5*(nrg-p->g)+rad);
    }
  if ( p->mu > .5*(p->g+nrg)-rad )
    {
      buff += fabs(rad-.5*(nrg-p->g));
    }
#else /*S=1/2*/
  double rad = sqrt( .25*(p->g-nrg)*(p->g-nrg)+.5625*(p->J)*(p->J)*(p->V)*(p->V) );
  if ( p->mu > .5*(p->g+nrg)+rad )
    {
      buff += (rad-.5*(p->g-nrg));
    }
  if ( p->mu > .5*(p->g+nrg)-rad )
    {
      buff += (rad+.5*(p->g-nrg));
    }
#endif /*USECASEJ32*/
  return buff/2./rad;
}


double rho_m_k(double x,       
	       double y,
	       double z,
	       void * params )
{
  struct Parameters * p = params;
  double nrg = EN(x,y,z);
#ifdef USECASEJ32
  double rad = sqrt( .25*(p->g-nrg)*(p->g-nrg)+14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
#else /*S=1/2*/
  double rad = sqrt( .25*(p->g-nrg)*(p->g-nrg)+.5625*(p->J)*(p->J)*(p->V)*(p->V) );
#endif /*USECASEJ32*/
  double buff = 0.;
  if ( p->mu > .5*(p->g+nrg)+rad )
    {
      buff += fabs(.5*(p->g-nrg)+rad);
    }
  if ( p->mu > .5*(p->g+nrg)-rad )
    {
      buff += fabs(rad-.5*(p->g-nrg));
    }
  return buff/2./rad;
}


/* this was experimental and appears to be incorrect */
double
rho_hyb( double x,      /*S12*/
	 void * param)
{
  struct Parameters * p = param;
#ifdef AFFINE
  return .5*(p->V) * dirac(x-(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g)) * fabs( (1-(x-p->g)*(x-p->g)/(p->J)/(p->J)/(p->V)/(p->V))/(x-p->g) );
#else
  return .375*(p->V) * dirac(x-.5625*(p->J)*(p->J)*(p->V)*(p->V)/(x-p->g)) * fabs( (1-(x-p->g)*(x-p->g)/.5625/(p->J)/(p->J)/(p->V)/(p->V))/(x-p->g) );
#endif
}





/********************************* Eigenval x occupation functions ******************************************/





double
eigen_occup( double x,      
	     void * param )
{
  struct Parameters * p = param;
  double ave = .5*(p->g+x) ;
  double buff = 0.;
#ifdef USECASEJ32
  double rad = sqrt( .25*(p->g-x)*(p->g-x)+14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
#else /*S=1/2*/
  double rad = sqrt( .25*(x-p->g)*(x-p->g) + .5625*(p->J)*(p->J)*(p->V)*(p->V) );
#endif /*USECASEJ32*/
  if (p->mu > ave-rad)
    buff = (ave-rad)*dirac(x) ;
  if (p->mu > ave+rad)
    buff += (ave+rad)*dirac(x) ;
  return buff;
}


double
eigen_occup_k( double x,      
	       double y,
	       double z,
	       void * param )
{
  struct Parameters * p = param;
  double nk = EN(x,y,z);
  double ave =.5*(p->g+nk) ;
  double buff = 0.;
#ifdef USECASEJ32
  double rad = sqrt( .25*(p->g-nk)*(p->g-nk)+14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
#else /*S=1/2*/
  double rad = sqrt( .25*(nk-p->g)*(nk-p->g) + .5625*(p->J)*(p->J)*(p->V)*(p->V) );
#endif /*USECASEJ32*/
  if (p->mu > ave-rad)
    buff = ave-rad ;
  if (p->mu > ave+rad) /* removed jun12 to check eigenocc calculation, but this affects only areas only where nc->filled. Otherwise, incorrect  */
    buff += ave+rad ;
  return buff;
}





/**************************************** SCE's and SCE-callers *******************************************/




#ifdef USECASEJ32


double
sceV(double x,      /*J32*/
     void * param)
{
  struct Parameters * p = param;
  int step;
  double delta = sqrt(.25*(x-p->g)*(x-p->g)+14.0625*p->J*p->J*p->V*p->V);
  double gpm = p->g/(-2.) + p->mu;
  if ( x/2. > gpm-delta )
    step=1;
  else
    step=0;
  if ( x/2. > gpm+delta )
    step--;
  if ( step < 0 )
    {
      return -1.875*dirac(x)/delta; 
    }
  if ( step > 0 )
    {
      return 1.875*dirac(x)/delta; 
    }  
  return 0.;
}

double
sce_k(double x,      /*J32*/
      double y,
      double z,
      void * param)
{
  struct Parameters * p = param;
  int step;
  double nrg = -2.*(cos(x)+cos(y)+cos(z));
  double gpm = p->g/(-2.) + p->mu;
  double delta = sqrt( .25*(nrg-p->g)*(nrg-p->g)+14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  if ( nrg/2. > gpm-delta )
    step=1;
  else
    step=0;
  if ( nrg/2. > gpm+delta )
    step--;
  if ( step < 0 )
    {
      return -1.875/delta;
    }
  if ( step > 0 )
    {
      return 1.875/delta; 
    }  
  return 0.;
}



#else /* S12 */


double
sceV(double x,      /*S12*/
     void * param)
{
  struct Parameters * p = param;
  int step;
  /* I will be using only p->J under the assumption that elsewhere I've set p->J=JV,p->V=1
   * I haven't accounted for half-filling here because I might want to look at other cases */
  double delta = sqrt( .25*(x-p->g)*(x-p->g) + .5625*(p->J)*(p->J) );
  double gpm = p->g/(-2.) + p->mu;
  if ( x/2. > gpm-delta )
    step=1;
  else
    step=0;
  if ( x/2. > gpm+delta )
    step--;
  if ( step < 0 )
    {
      //return -.375*rho_c(x,param)/delta;
      //return -.375*rho_m(x,param)/delta; /* returns closest value to sce_k */
      return -.375*dirac(x)/delta;
    }
  if ( step > 0 )
    {
      //return .375*rho_c(x,param)/delta;
      //return .375*rho_c(x,param)/delta;  /* returns closest value to sce_k */
      return .375*dirac(x)/delta;
    }  
  return 0.;
}

double
sce_k(double x,      /*S12*/
      double y,
      double z,
      void * param)
{
  struct Parameters * p = param;
  int step;
  /* I will be using only p->J under the assumption that elsewhere I've set p->J=JV,p->V=1
   * I haven't accounted for half-filling here because I might want to look at other cases */
  double nrg = -2.*(cos(x)+cos(y)+cos(z));
  double gpm = p->g/(-2.) + p->mu;
  double delta = sqrt( .25*(nrg-p->g)*(nrg-p->g) + .5625*(p->J)*(p->J) );
  if ( nrg/2. > gpm-delta )
    step=1;
  else
    step=0;
  if ( nrg/2. > gpm+delta )
    step--;
  if ( step < 0 )
    {
      return -.375/delta;
    }
  if ( step > 0 )
    {
      return .375/delta;
    }
  return 0.;
}



#endif /*USECASEJ32*/




/***********************************  msplit & msolve *****************************************/




int
msplit(struct Parameters * par,
       double lb,
       double ub,
       double * lmu,
       double * umu)
{
  double inc;
  double test;
  char * code;

  /* gamma = 0 is a special case */

  if ( NEARZERO(par->g) )
    {
#ifdef DEBUG
      fprintf(par->fout,"    : lambda is near zero\n");
#endif
      *lmu = msolve(par,lb,0.,0.);
      *umu = *lmu;
      return 0;
    }

  /* gamma != 0 */

  test = quad(rho_m,par,LOWERINT,par->g,0.,REL_ERROR,&code) - NIMP;
#ifdef DEBUG
  fprintf(par->fout,"    : %1.5e ",test);
#endif
  if ( code!=NULL )
    {
      free(code);
      code=NULL;
      return -1;
    }
  if ( NEAR(test,0.,PREC/10.) )
    {
      /* We've found the plateau                         *
       * split the search domain and call msolve         */
      if ( test<0. )
	{
#ifdef DEBUG
	  fprintf(par->fout,"plateau offsets +1,0 ");
#endif
	  *lmu = msolve( par, lb, par->g, 1. ); /* add offset */
	  *umu = msolve( par, par->g, ub, 0. ); /* no offset */
	}
      else
	{
#ifdef DEBUG
	  fprintf(par->fout,"plateau offsets 0,-1 ");
#endif
	  *lmu = msolve( par, lb, par->g, 0. ); /* no offset */
	  *umu = msolve( par, par->g, ub, -1. ); /* subtract offset */
	}
    }
  else
    {
      if ( test > 0. )
	{
#ifdef DEBUG
      fprintf(par->fout,"lower interval ");
#endif
	  *lmu = msolve( par, lb, par->g, 0. );
	  *umu = *lmu;
	}
      else
	{
#ifdef DEBUG
      fprintf(par->fout,"upper interval ");
#endif
	  *lmu = msolve( par, par->g, ub, 0. );
	  *umu = *lmu;	  
	}
    }
      
  return 0;

}


double
msolve( struct Parameters * par,
	const double lb,
	const double ub,
	const double offset )
{
  double inc;
  double mvalues[3];
  double resultsMID, resultsHI;
  double eta = PREC/10.;
  char * code;

  mvalues[LOW] = lb;
  mvalues[HI] = ub;
  mvalues[MID] = (lb+ub)/2.;
  
  resultsHI = quad(rho_m,par,LOWERINT,mvalues[HI],0.,REL_ERROR,&code) - NIMP + offset*eta;
#ifdef DEBUG
  fprintf(par->fout,": %1.5e ",resultsHI);
#endif
  
  if ( code!=NULL )
    {
      free(code);
      code=NULL;
      return 2.*UPPER;
      }

  if ( resultsHI < 0. )
    {
#ifdef DEBUG
      fprintf(par->fout,"fail\n");
#endif
      return 2*UPPER;
    }
#ifdef DEBUG
  else
    {
      fprintf(par->fout,"\n");
    }
#endif

  inc = ub - lb;

  while ( inc > PREC )
    {

      resultsMID = quad(rho_m,par,LOWERINT,mvalues[MID],0.,REL_ERROR,&code) - NIMP + offset*eta;
#ifdef DEBUG
      fprintf(par->fout,"    :       : [%1.5e,%1.5e] (%1.5e) %1.5e\n",mvalues[LOW],mvalues[HI],mvalues[MID],resultsMID);
#endif
      
      if ( code!=NULL )
	{
          free(code);
	  code=NULL;
	  }

      /* if ( SIGNCHANGE(results[HI],results[MID]) ) */
      if ( SIGNCHANGE(resultsHI,resultsMID) )
	{
	  mvalues[LOW] = mvalues[MID];
	  mvalues[MID] = (mvalues[MID]+mvalues[HI])/2.;
	}
      else
	{
	  /* results[HI] = results[MID]; */
	  resultsHI = resultsMID;
	  mvalues[HI] = mvalues[MID];
	  mvalues[MID] = (mvalues[LOW]+mvalues[MID])/2.;
	}
      //inc /= 2.;
      inc = mvalues[HI]-mvalues[LOW];
    } 

  return (mvalues[LOW]+mvalues[HI])/2.;
  
}




int
msplit_k(struct Parameters * par,
	 double lb,
	 double ub,
	 double * lmu,
	 double * umu)
{
  double inc;
  double test;
  char * code;

  /* gamma = 0 is a special case */

  if ( NEARZERO(par->g) )
    {
#ifdef DEBUG
      fprintf(par->fout,"    : lambda is near zero\n");
#endif
      *lmu = msolve_k(par,lb,0.,0.);
      *umu = *lmu;
      return 0;
    }

  /* gamma != 0 */

  par->mu = par->g;
  test = SUMMOR(rho_m_k,par,SUMPREC) - NIMP;
#ifdef DEBUG
  fprintf(par->fout,"    : %1.5e ",test);
#endif
  if ( NEAR(test,0.,PREC/10.) )
    {
      /* We've found the plateau                         *
       * split the search domain and call msolve         */
      if ( test<0. )
	{
#ifdef DEBUG
	  fprintf(par->fout,"plateau offsets +1,0 ");
#endif
	  *lmu = msolve_k( par, lb, par->g, 1. ); /* add offset */
	  *umu = msolve_k( par, par->g, ub, 0. ); /* no offset */
	}
      else
	{
#ifdef DEBUG
	  fprintf(par->fout,"plateau offsets 0,-1 ");
#endif
	  *lmu = msolve_k( par, lb, par->g, 0. ); /* no offset */
	  *umu = msolve_k( par, par->g, ub, -1. ); /* subtract offset */
	}
    }
  else
    {
      if ( test > 0. )
	{
#ifdef DEBUG
      fprintf(par->fout,"lower interval ");
#endif
	  *lmu = msolve_k( par, lb, par->g, 0. );
	  *umu = *lmu;
	}
      else
	{
#ifdef DEBUG
      fprintf(par->fout,"upper interval ");
#endif
	  *lmu = msolve_k( par, par->g, ub, 0. );
	  *umu = *lmu;	  
	}
    }
      
  return 0;

}


double
msolve_k( struct Parameters * par,
	  const double lb,
	  const double ub,
	  const double offset )
{
  double inc;
  double mvalues[3];
  double resultsMID, resultsHI;
  double eta = PREC/10.;
  char * code;

  mvalues[LOW] = lb;
  mvalues[HI] = ub;
  mvalues[MID] = (lb+ub)/2.;
  
  par->mu = mvalues[HI];
  resultsHI = SUMMOR(rho_m_k,par,SUMPREC) - NIMP + offset*eta;
#ifdef DEBUG
  fprintf(par->fout,": %1.5e ",resultsHI);
#endif

  if ( resultsHI < 0. )
    {
#ifdef DEBUG
      fprintf(par->fout,"fail\n");
#endif
      return 2*UPPER;
    }
#ifdef DEBUG
  else
    {
      fprintf(par->fout,"\n");
    }
#endif

  inc = ub - lb;

  while ( inc > PREC )
    {

      par->mu = mvalues[MID];
      resultsMID = SUMMOR(rho_m_k,par,SUMPREC) - NIMP + offset*eta;
#ifdef DEBUG
      fprintf(par->fout,"    :       : [%1.5e,%1.5e] (%1.5e) %1.5e\n",mvalues[LOW],mvalues[HI],mvalues[MID],resultsMID);
#endif

      if ( SIGNCHANGE(resultsHI,resultsMID) )
	{
	  mvalues[LOW] = mvalues[MID];
	  mvalues[MID] = (mvalues[MID]+mvalues[HI])/2.;
	}
      else
	{
	  resultsHI = resultsMID;
	  mvalues[HI] = mvalues[MID];
	  mvalues[MID] = (mvalues[LOW]+mvalues[MID])/2.;
	}
      inc = mvalues[HI]-mvalues[LOW];
    } 

  return (mvalues[LOW]+mvalues[HI])/2.;
  
}




int
mgsolve( struct Parameters * par,
	 double lb,
	 double ub )
{
  double inc;
  double mvalues[3];
  double results[3];
  double buffer;
  char * code;
  
  mvalues[LOW] = lb;
  mvalues[HI] = ub;
  mvalues[MID] = (lb+ub)/2.;

#  ifdef DEBUG
  fprintf(par->fout,"[%1.3e,%1.3e] ",mvalues[LOW],mvalues[HI]);
#  endif

  par->g = mvalues[HI];
  par->mu=par->g;
#  ifdef DEBUG
  buffer = quad(rho_m,par,lb,par->mu,0.,REL_ERROR,&code);
  results[HI] = buffer - NIMP;
  fprintf(par->fout,"%1.6e - %f = %1.6e\n",buffer,NIMP,results[HI]);
#  else
  results[HI] = quad(rho_m,par,lb,par->mu,0.,REL_ERROR,&code) - NIMP;
#  endif
  
  if ( code!=NULL )
    {
#  ifdef DEBUG
      fprintf(par->fout,"code: %s\n",code);
#  endif
      free(code);
      code=NULL;
    }
  if ( results[HI] > 0 ) /***** was > ******/
    {
#  ifdef DEBUG
      fprintf(par->fout,"NIMP not found\n");
#  endif
      return -1;
    }

  inc = ub - lb;

  while ( inc > PREC )
    {

      par->g = mvalues[MID];
      par->mu=par->g;
      results[MID] = quad(rho_m,par,lb,par->mu,0.,REL_ERROR,&code) - NIMP;
#  ifdef DEBUG
      buffer = quad(rho_m,par,lb,par->mu,0.,REL_ERROR,&code);
      results[MID] = buffer - NIMP;
      fprintf(par->fout,"[%1.3e,%1.3e] mid: %1.6e, %1.6e - %f = %1.6e\n",mvalues[LOW],mvalues[HI],mvalues[MID],buffer,NIMP,results[MID]);
#  else
      results[MID] = quad(rho_m,par,lb,par->mu,0.,REL_ERROR,&code) - NIMP;
#  endif
      if ( code!=NULL )
	{
#  ifdef DEBUG
	  fprintf(par->fout,"code: %s\n",code);
#  endif
	  free(code);
	  code=NULL;
	}

      if ( SIGNCHANGE(results[HI],results[MID]) )
	{
	  results[LOW] = results[MID];
	  mvalues[LOW] = mvalues[MID];
	  mvalues[MID] = (mvalues[MID]+mvalues[HI])/2.;
	}
      else
	{
	  results[HI] = results[MID];
	  mvalues[HI] = mvalues[MID];
	  mvalues[MID] = (mvalues[LOW]+mvalues[MID])/2.;
	}
      /*inc /= 2.;*/
      inc = results[HI]-results[LOW];
    } 
  
  par->g = (mvalues[LOW]+mvalues[HI])/2.;
  par->mu = par->g;
  return 0;
}
