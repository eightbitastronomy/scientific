/**********************************************************************************
 *  bulkcore.c - Library for phase bulk/nonboundary calculations
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
 *******************************************************************************************************************************
 *
 * 2018 aug 22: mu_finder and glam_finder now act on copies of parameters; arguments are the same. 
 *              Added exit criterion to mu_finder and glam_finder.
 * 2018 aug 20: attempts to correct Kondo energy have failed; implementing free_energy_kondo2 as a remedy.
 * 2018 aug 16: debug3 added to free_energy_afm for printing
 * 2018 aug 13: 
 *   o  kondo free_energy altered to calculate k-summation terms only
 *   o  afm free energy recoded. New fctn free_energy_afm handles overall terms, free_energy_afm_band_wrapper handles
 *      conduction-band summation & integration.
 *
 *
 **************************************************************************************************************************/

#include "bulkcore.h"
#include "string.h"





extern double SUMMOR(double (*f)(double,double,double,void *),void *,double);
extern double quad(double (*f)(double,void *),void * p,double lb,double ub,double abs,double rel,char ** ecode);







/**************************  nc **************************************/



double occupation_kondo_c( double x,
			   double y,
			   double z,
			   void * param )
{
  struct Parameters * p = param;
  double nk = EN(x,y,z);
  double emg = ( nk-p->g ) / 2. ;
  double rad = sqrt(  .25*(p->g-nk)*(p->g-nk) + 14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  /* multiplied by 4 to account for all 4 spins */
  return 2.*( fabs(rad+emg)*FERMI(rad+.5*(p->g+nk),p->mu,KB,p->T) + fabs(rad-emg)*FERMI(.5*(p->g+nk)-rad,p->mu,KB,p->T) ) / rad ;
}




double occupation_afm_c_wrapper( double x,
				 void * param )
{
  struct Pwrapper pw;
  pw.e = x;
  pw.p = param;
  return FERMI( x, (pw.p)->mu, KB, (pw.p)->T )  *  SUMMOR( heatkernel_reduced,&pw,SUMPREC ) ;
}





/************************** nf ******************************************/




double occupation_kondo_f( double x,
			   double y,
			   double z,
			   void * param )
{
  struct Parameters * p = param;
  double nk = EN(x,y,z);
  double gme = ( p->g-nk ) / 2. ;
  double rad = sqrt(  gme*gme + 14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  /* multiplied by 4 to account for all 4 spins */
  return 2.*( fabs(rad+gme)*FERMI(rad+.5*(p->g+nk),p->mu,KB,p->T) + fabs(rad-gme)*FERMI(.5*(p->g+nk)-rad,p->mu,KB,p->T) ) / rad ;
}






double occupation_afm_f( double g,
			 double mu,
			 double j,
			 double s,
			 double t )
{
  return FERMI(g+1.5*j*s,mu,KB,t)
    + FERMI(g-1.5*j*s,mu,KB,t)
    + FERMI(g+.5*j*s,mu,KB,t)
    + FERMI(g-.5*j*s,mu,KB,t);
}





/************************** lindhard ********************************/




double sce_kondo( double x,
		  double y,
		  double z,
		  void * param )
{
  struct Parameters * p = param;
  double nk = EN(x,y,z);
  double gme = ( p->g-nk ) / 2. ;
  double rad = sqrt(  gme*gme + 14.0625*(p->J)*(p->J)*(p->V)*(p->V) );
  return 1.875 * ( FERMI(.5*(p->g+nk)-rad,p->mu,KB,p->T)-FERMI(.5*(p->g+nk)+rad,p->mu,KB,p->T) ) / rad ;
}



double sce_afm_m( double g,
		  double mu,
		  double j,
		  double s,
		  double t )
{
  return 1.5 * ( FERMI(g-1.5*j*s,mu,KB,t)-FERMI(g+1.5*j*s,mu,KB,t) )
    + .5 * ( FERMI(g-.5*j*s,mu,KB,t)-FERMI(g+.5*j*s,mu,KB,t) )       ;
}



double sce_afm_s_wrapper( double x,
			  void * param )
{
  struct Pwrapper pw;
  pw.e = x;
  pw.p = param;
  return FERMI( x, (pw.p)->mu, KB, (pw.p)->T )  *  SUMMOR( sce_afm_s,&pw,SUMPREC ) ;
}




double sce_afm_s( double x,
		  double y,
		  double z,
		  void * param )
{
  double e = ( (struct Pwrapper *)param )->e ;
  struct Parameters * p = ( (struct Pwrapper *)param )->p ;
  /*double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double enk = e - nk;
  double enkq = e - nkq;
  double enk2q = e - nk2q;
  double enk3q = e - nk3q;*/
  double enk = e - EN(x,y,z);
  double enkq = e - EN(x+p->qx,y+p->qy,z+p->qz);
  double enk2q = e - EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double enk3q = e - EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);
  double B = 1.5*enkq*enk + 2.*enk*enk3q + 1.5*enk2q*enk3q - 2.25*jmsq ;
  double A = enk * enkq * enk2q * enk3q
    - jmsq * ( .75*enk*enkq + enk*enk3q + .75*enk2q*enk3q )
    + .5625 * jmsq*jmsq;
  double dAdw = (enk+enk3q)*(enkq*enk2q-1.75*jmsq) + (enkq+enk2q)*(enk*enk3q-.75*jmsq) ;
  return -1. * p->J * p->magM * B / dAdw * exp( A*A/dAdw/dAdw/(-2.)/(p->eta)/(p->eta) ) / sqrt(2.*M_PI) / (p->eta) ;
}





/************************** heat kernel *****************************/




double
heatkernel ( double x,
	     double y,
	     double z,
	     void * param )
{
  /* energy will enter through a struct { double, struct Parameters * } */
  double e = ( (struct Pwrapper *)param )->e ;
  struct Parameters * p = ( (struct Pwrapper *)param )->p ;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double A = ( e - nk ) * ( e - nkq ) * ( e - nk2q ) * ( e - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM) * ( .75*(e-nk)*(e-nkq) + (e-nk)*(e-nk3q) + .75*(e-nk2q)*(e-nk3q) )
    + .5625 * (p->J)*(p->J)*(p->J)*(p->J)*(p->magM)*(p->magM)*(p->magM)*(p->magM);
  double dAdw = ( e - nk ) * ( e - nkq ) * ( e - nk2q )
    + ( e - nk ) * ( e - nkq ) * ( e - nk3q )
    + ( e - nk ) * ( e - nk2q ) * ( e - nk3q )
    + ( e - nkq ) * ( e - nk2q ) * ( e - nk3q )
    - (p->J)*(p->J)*(p->magM)*(p->magM)*( 1.75*(2.*e-nk-nk3q) + .75*(2.*e-nkq-nk2q) );
  return exp( A*A/dAdw/dAdw/(-2.)/(p->eta)/(p->eta) ) / sqrt(2.*M_PI) / (p->eta) ;
}



double
heatkernel_reduced ( double x,
		     double y,
		     double z,
		     void * param )
{
  /* energy will enter through a struct { double, struct Parameters * } */
  double e = ( (struct Pwrapper *)param )->e ;
  struct Parameters * p = ( (struct Pwrapper *)param )->p ;
  double nk = EN(x,y,z);
  double nkq = EN(x+p->qx,y+p->qy,z+p->qz);
  double nk2q = EN(x+2.*p->qx,y+2.*p->qy,z+2.*p->qz);
  double nk3q = EN(x+3.*p->qx,y+3.*p->qy,z+3.*p->qz);
  double enk = e - nk;
  double enkq = e - nkq;
  double enk2q = e - nk2q;
  double enk3q = e - nk3q;
  double jmsq = (p->J)*(p->J)*(p->magM)*(p->magM);
  double A = enk * enkq * enk2q * enk3q
    - jmsq * ( .75*enk*enkq + enk*enk3q + .75*enk2q*enk3q )
    + .5625 * jmsq*jmsq;
  double dAdw = (enk+enk3q)*(enkq*enk2q-1.75*jmsq) + (enkq+enk2q)*(enk*enk3q-.75*jmsq) ;
  return exp( A*A/dAdw/dAdw/(-2.)/(p->eta)/(p->eta) ) / sqrt(2.*M_PI) / (p->eta) ;
}





/***************************** Free energy ***********************************************/


int mastercounter=0;


/* constant values in F must be added after this function returns */
double
free_energy_kondo ( double x,
		    double y,
		    double z,
		    void * param )
{
  struct Parameters * p = param;
  double nrg = EN(x,y,z);
  double epg = .5*(nrg+p->g);
  double rad = sqrt( .25*(p->g-nrg)*(p->g-nrg)+14.0625*p->J*p->J*p->V*p->V );
  double fplus = FERMI(epg+rad,p->mu,KB,p->T);
  double fminus = FERMI(epg-rad,p->mu,KB,p->T);
    
#ifdef DEBUG2
  double xfxp,fxlp=0.,mfxlp=0.,xfxm=0.,fxlm=0.,mfxlm=0.;
  xfxp = (epg+rad)*fplus;
  xfxm = (epg-rad)*fminus;
  if ( BETWEEN(fplus) )
    {
      fxlp = fplus*log(fplus);
    }
  if ( BETWEEN(fminus) )
    {
      fxlm = fminus*log(fminus);
    }
  if ( (1.-fplus) > 1.e-14 )
    {
      mfxlp = (1.-fplus)*log(1.-fplus);
    }
  if ( (1.-fminus) > 1.e-14 )
    {
      mfxlm = (1.-fminus)*log(1.-fminus);
    }
  if ( mastercounter == 1000000 )
    mastercounter = 0;
  if ( mastercounter++ % 100 == 0 )
    {
      fprintf(p->fout,"          (%1.4e,%1.4e): xfxp %1.4e %1.4e fxl %1.4e %1.4e mfxl %1.4e %1.4e\n",fplus,fminus,xfxp,xfxm,fxlp,fxlm,mfxlp,mfxlm);
      fflush(p->fout);
    }
  return xfxp + xfxm + KB*(p->T)*(fxlp+fxlm+mfxlp+mfxlm) ;
#else
  double fsum, flogsum;
  fsum = (epg+rad)*fplus;
  flogsum=0.;
  if ( NEARZERO(fplus)==0 )
    {
      flogsum += fplus*log(fplus);
    }
  fsum += (epg-rad)*fminus;
  if ( NEARZERO(fminus)==0 )
    {
      flogsum += fminus*log(fminus);
    }
  if ( NEARZERO(1.-fplus)==0 )
    {
      flogsum += (1.-fplus)*log(1.-fplus);
    }
  if ( NEARZERO(1.-fminus)==0 )
    {
      flogsum += (1.-fminus)*log(1.-fminus);
    }
  return fsum + KB*(p->T)*flogsum ;
#endif
}




double
free_energy_kondo2 ( double x,
		     double y,
		     double z,
		     void * param )
{
  struct Parameters * p = param;
  double nrg = EN(x,y,z);
  double epg = .5*(nrg+p->g);
  double rad = sqrt( .25*(p->g-nrg)*(p->g-nrg)+14.0625*p->J*p->J*p->V*p->V );
  double fsum, flogsum;

  //flogsum = log( 1. + exp(-1.*(epg+rad-p->mu)/KB/(p->T)) );
  flogsum = log1p(exp(-1.*(epg+rad-p->mu)/KB/(p->T)));
  //flogsum += log( 1. + exp(-1.*(epg-rad-p->mu)/KB/(p->T)) );
  flogsum += log1p(exp(-1.*(epg-rad-p->mu)/KB/(p->T)));
  
  return -1.*KB*(p->T)*flogsum ;
}



/* constant values in F must be added after this function returns */
double
free_energy_afm_band_wrapper ( double x,
			       void * param )
{
  struct Pwrapper pw;
  struct Parameters * r = param;
  pw.e = x;
  pw.p = param;
  double fx = FERMI(x,r->mu,KB,r->T);
  double total = 0.;
  if ( NEARZERO(fx) )
    return SUMMOR( heatkernel_reduced, &pw, SUMPREC ) * KB * r->T * (1.-fx)*log(1.-fx);
  if ( NEARZERO(1.-fx) )
    return SUMMOR( heatkernel_reduced, &pw, SUMPREC ) * ( x*fx+ KB*(r->T)*fx*log(fx) );
  return SUMMOR( heatkernel_reduced, &pw, SUMPREC ) * (  x*fx
							 + KB*(r->T)*( fx*log(fx)
								       + (1.-fx)*log(1.-fx) )
							 );
}



double
free_energy_afm(void * param)
{
  struct Parameters * p = param;
  double cbands;
  double constants;
  double fbands = 0.;
  double feigs[4];
  double ffermis[4];
  int i;
  char * ecode = NULL;
#ifdef DEBUG3
  fprintf(p->fout,"        [     entering AFM calc\n");
  fflush(p->fout);
#endif
  cbands = quad( free_energy_afm_band_wrapper,p,LOWER,UPPER,0,REL_ERROR,&ecode );
  if (ecode)
    {
      free(ecode);
      ecode=NULL;
    }
#ifdef DEBUG3
  fprintf(p->fout,"        [     Cbands %1.4e\n",cbands);
  fflush(p->fout);
#endif
  constants = (p->J)*(p->magS)*sce_afm_m(p->g,p->mu,p->J,p->magS,p->T);
#ifdef DEBUG3
  fprintf(p->fout,"        [     Constants %1.4e\n",constants);
  fflush(p->fout);
#endif
  feigs[0] = p->g + 1.5*(p->J)*(p->magS);
  feigs[1] = p->g - 1.5*(p->J)*(p->magS);
  feigs[2] = p->g + .5*(p->J)*(p->magS);
  feigs[3] = p->g - .5*(p->J)*(p->magS);
  ffermis[0] = FERMI( feigs[0],p->mu,KB,p->T );
  ffermis[1] = FERMI( feigs[1],p->mu,KB,p->T );
  ffermis[2] = FERMI( feigs[2],p->mu,KB,p->T );
  ffermis[3] = FERMI( feigs[3],p->mu,KB,p->T );
  for ( i=0; i<4; i++ )
    {
#ifdef DEBUG3
      fprintf(p->fout,"        [     %d feig %1.4e, ffermi %1.4e\n",i,feigs[i],ffermis[i]);
      fflush(p->fout);
#endif
      fbands += feigs[i]*ffermis[i];
      if ( NEARZERO(ffermis[i])==0 )
	{
#ifdef DEBUG3
	  fprintf(p->fout,"        [     %d log ffermi %1.4e, kT %1.4e\n",log(ffermis[i]), KB*p->T);
	  fflush(p->fout);
#endif
	  fbands += KB*(p->T)*ffermis[i]*log(ffermis[i]);
	}
      if ( NEARZERO(1.-ffermis[i])==0 )
	{
#ifdef DEBUG3
	  fprintf(p->fout,"        [     log(1-ffermi) %1.4e, kT %1.4e\n",log(1.-ffermis[i]), KB*p->T);
	  fflush(p->fout);
#endif
	  fbands += KB*(p->T)*(1.-ffermis[i])*log(1.-ffermis[i]);
	}
    }
  return cbands + constants - p->g + fbands;
}





/*******************   mu finder for bulk   ***************************/





double
mu_finder( double (*mf)( void *, double ),
	   void * par,
	   double targetnc )
{
  struct Parameters p;
  char * ecode = NULL;
  double lower = LOWER;
  double target;
  double res;
  double upper = UPPER;

  memcpy(&p,par,sizeof(struct Parameters));
  
#ifdef DEBUG
  //  fprintf(ip->fout,"------->>>> [%1.4e,%1.4e] [%1.4e,%1.4e] nc *%f* \n",lower,ip->mu,ip->mu,upper,ip->nc);
  //fflush(ip->fout);
#endif  
  do
    {
      res = mf( &p,targetnc );
#ifdef DEBUG
      //fprintf(ip->fout,"------->>>> [%1.4e,%1.4e] [%1.4e,%1.4e] nc %1.4e\n",lower,ip->mu,ip->mu,upper,res);
      //fflush(ip->fout);
#endif
      if (NEAR(res,targetnc,PREC))
	return p.mu;
      if (res < targetnc)
	{
	  lower=p.mu;
	}
      else
	{
	  upper=p.mu;
	}      
      p.mu = (lower+upper)/2.;

    } while (NEAR(upper,lower,PREC/2.)==0);

  return p.mu;
}




double mu_wrapper_afm( void * par,
		       double trn  )
{
  struct Parameters * p = par;
  char * ecode = NULL;
  double result = quad( occupation_afm_c_wrapper,par,LOWER,UPPER,0,REL_ERROR,&ecode );
  if (ecode)
    {
      free(ecode);
    }
  p->nc = result;
  return result;
}



double mu_wrapper_kondo( void * par,
			 double trn  )
{
  return SUMMOR( occupation_kondo_c,par,SUMPREC );
}





double
glam_finder( double (*gf)(void *, double),
	     void * par,
	     double nimp )
{
  struct Parameters p ;
  double lower = LOWER;
  double target;
  double res;
  double upper = UPPER;

  memcpy(&p,par,sizeof(struct Parameters));

#ifdef DEBUG
  //fprintf(p->fout,"-------^^^^^ [%1.4e,%1.4e][%1.4e,%1.4e] ==> .25\n",lower,target,target,upper);
  //fflush(p->fout);
#endif
  target = lower;
  res = gf( &p,target );
  if ( res < nimp )
    return target;

  do
    {
      target = (lower+upper)/2.;
      res = gf( &p,target );
#ifdef DEBUG
      //fprintf(p->fout,"-------^^^^^ [%1.4e,%1.4e][%1.4e,%1.4e] ==> %1.4e\n",lower,target,target,upper,res);
      //fflush(p->fout);
#endif
      if ( NEAR(res,nimp,PREC) )
	return target;
      if (res > nimp)
	{
	  lower=target;
	}
      else
	{
	  upper=target;
	}

    } while (NEAR(upper,lower,PREC/2.)==0);
  return target;
}



double glam_wrapper_afm( void * par,
			 double trg )
{
  struct Parameters * p = par;
  p->g = trg;
  return occupation_afm_f( p->g,p->mu,p->J,p->magS,p->T )  ;
}




double glam_wrapper_kondo( void * par,
			   double trg )
{
  struct Parameters * p = par;
  p->g = trg;
  return SUMMOR( occupation_kondo_f,p,SUMPREC );
}




/**********************  minimization / filtration ****************************/






int
filter_minimum( double * record,
		int recordsz,
		int fields,
		int key             )
{
  int i;
  int best;
  double bestval = 0.;

  for (i=0;i<recordsz;i++)
    {
      if ( record[fields*i+key]<bestval )
	{
	  bestval = record[fields*i+key] ;
	  best = i;
	}
    }

  return best;
}
