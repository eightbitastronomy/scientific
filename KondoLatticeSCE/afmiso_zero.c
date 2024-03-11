/**********************************************************************************
 *  afmiso_c - T=0 calculator for Isotropic-model SDW energy
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
 ***********************************************************************************
 *
 *   For grouping dispatches by nc, iteration in q
 *   Updates:
 *     29 jun 2018,
 *        o  domain on integration modified to +/-6t +/- JM or +/-10. 
 *        o  mu_finder & mu_stepper modified to accept dynamic bounds
 *        o  mu_finder altered to short-circuit if test-result is within PREC/10 of target nc
 *        o  paramagnetic energy calcs produce zeroes until I can fix the bare.c functions
 *        o  removed any traces of krandom library
 *     6 jun 2018, 
 *        o  removed incommensurate q's. krand library no longer needed.
 *        o  combined "lump" version of code with this one. Use -DLUMP in makefile if desired.
 *        o  reinstated J-iteration functionality.
 *     25 may 2018, COMMENSURATE_Q flag was misapplied. Corrected.
 *     23 may 2018,
 *        o  inclusion of incommensurate q values 
 *        o  changed mesh data directory to be relative. Uncommented stop_ninterp at end.
 *     2 may 2018, repurposed from sdw_zero_general.c. Iteration is still allowed over J, but is not being used,
 *       and printout of J has been removed.
 *       mu_finder has fixed upper/lower bounds now, since JM is currently fixed.
 *
 *************************************************************************************************************************/



/*************************************************
*  options to implement, if possible
*     --nlow, nhi, ninc
*     --jlow, jhi, jinc
*     --qsplit vs --q
*     --qcom vs --qincom vs --qall
*     --sumprec, --prec, --relerr
*     --nint
*     --eta
*     --simp, --trap
*     --J
*     --lump vs --loopnj vs --loopq
*************************************************/


#include <stdio.h>
#include <stdarg.h>
#include "nint.h"
#include "sm.h"
#include "fio.h"
#include "bare.h"
#include "afmisozerocore.h"


const char dirname[8] = "./input\0";
const char meshname[17] = "dos_mesh_t_1.txt\0";


#define LOWERJ             ( 1. )
#define COUNTJ             ( 1 )
#define INCJ               ( .5 )
#define LOWERNC            (3.5)
#define INCNC              ( .2 )
#define COUNTNC            ( 1 )
#define QSPLIT             ( 1 )
#define ETA                ( .01 )
#define PRIMARYFIELD       ( 6 )
#define SECONDARYFIELD     ( 3 )
#define COUNTFIELD         ( PRIMARYFIELD + SECONDARYFIELD )
#ifdef LUMP
#  define DATAFIELD        ( 5 )
#else
#  define DATAFIELD        ( 2 )
#endif


#define FINDMU             mu_stepper
  

/***  PRIMARY datafields: ***/
enum { MU,          /* chemical potential */
       N,           /* cond occupation */
       EP,          /* paramagnetic energy */
       EN,          /* cond occupation in PM phase */
       EM,          /* chemical potential in PM phase */
       R    };      /* SCE return value */

/*** SECONDARY datafields ***/
enum { S,           /* Magnetization "S" -- conduction */
       G,           /* Gamma or lambda -- langrange multiplier */
       E    };      /* Free energy */




void * masterRHS( void * par, void * in, void * out );
void egress(int ac, ... );





int main()
{

  int i,j,k,l,m,n;
  int ret;
  int qtotal = ((QSPLIT*QSPLIT*QSPLIT)+(3*QSPLIT*QSPLIT)+(2*QSPLIT))/6;
  int datacount;
  double incq = ((2*M_PI)/((double)QSPLIT))-PREC;
  double * dpts, * meshx, *meshy;
  int meshn;
  double * results = NULL;
  char * errorcode = NULL;
  FILE * input;
  struct Parameters p;
  struct d_Function df;


#ifdef DEBUG
  p.fout=stdout;
#endif


  /* ahead of calculations, print out some information */

  printf("#J=3/2 |M|=%f\n",MAGM);
  fprintf(stderr,"J=3/2 |M|=%f\n",MAGM);
  printf("#jlower=%f #NIMP is %f\n",JLOWER,NIMP);
  fprintf(stderr,"jlower=%f NIMP is %f\n",JLOWER,NIMP);
  printf("#Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stderr,"Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  printf("#%d divisions of q: %d calculations\n",QSPLIT,qtotal);
  fprintf(stderr,"%d divisions of q: %d calculations\n",QSPLIT,qtotal);

  
  /* open output files */

  printf("##x       y        z        nc   J            mu            nc2          EP            EPN          EPM          ret           S            g             E\n");


  /* nint preparation */

  meshn = get_mesh(dirname,meshname,&dpts);
  if ( meshn<1 )
    {
      fprintf(stderr,"Mesh problem: %d\n",meshn);
      if (dpts)
	free(dpts);
    }
  ret = mesh_copy( &meshx, &meshy, dpts, meshn );
  if ( ret < 0 )
    {
      if (dpts)
	free(dpts);
    }
  if (dpts)
    free(dpts);
  start_ninterp(meshn,meshx,meshy);
  meshx=NULL;
  meshy=NULL;
  start_nint(NINTINT);

  
  /* informational output */

  fprintf(stderr,"q-split will result in %d iterations...\n",qtotal);
#ifdef LUMP
  dpts = malloc(qtotal*COUNTNC*COUNTJ*DATAFIELD*sizeof(double));
  l=0;
  for (m=0;m<COUNTJ;m++)
    {
      for (n=0;n<COUNTNC;n++)
	{
	  for (i=0;i<QSPLIT;i++)
	    {
	      for (j=0;j<=i;j++)
		{
		  for (k=0;k<=j;k++)
		    {
		      dpts[l++]=LOWERJ+((double)m)*INCJ;
		      dpts[l++]=LOWERNC+((double)n)*INCNC;
		      dpts[l++]=((double)i)*incq-M_PI;
		      dpts[l++]=((double)j)*incq-M_PI;
		      dpts[l++]=((double)k)*incq-M_PI;
		    }
		}
	    }
	}
    }
  datacount = qtotal*COUNTNC*COUNTJ;
#else  
  dpts = malloc(COUNTNC*COUNTJ*DATAFIELD*sizeof(double));
  k = 0;
  for (i=0;i<COUNTJ;i++)
    {
      for (j=0;j<COUNTNC;j++)
	{
	  dpts[k++]=LOWERJ+((double)i)*INCJ;
	  dpts[k++]=LOWERNC+((double)j)*INCNC;
	}
    }
  datacount = COUNTNC*COUNTJ;
#endif
  
  fprintf(stderr,"%d(J) x %d(nc) x %d(q) = %d records\n",COUNTJ,COUNTNC,qtotal,COUNTJ*COUNTNC*qtotal);
  fprintf(stderr,"Looping...\n");

  p.magM = MAGM;
  p.eta = ETA;
  df.f   = masterRHS;
  df.param = &p;


  fflush(stderr);
  fflush(stdout);

  
  /* dispatch loop */

#ifndef LUMP
  for (i=0;i<QSPLIT;i++)
    {
      for (j=0;j<=i;j++)
	{
	  for (k=0;k<=j;k++)
	    {
	      p.qx=((double)i)*incq-M_PI;
	      p.qy=((double)j)*incq-M_PI;
	      p.qz=((double)k)*incq-M_PI;
#endif

	      /* dispatch */
	      
	      dispatchv2( df, datacount, dpts, DATAFIELD*sizeof(double), &results, COUNTFIELD*sizeof(double), egress );

	      /* dump results */
	      
	      for (l=0;l<datacount;l++)
		{
		  /*      qx    qy    qz    nc    J     mu    nc    ep    en    em    ret   S     g     E       */
		  printf("%1.5f %1.5f %1.5f %1.2f %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n",
#ifdef LUMP
			 dpts[l*DATAFIELD+2],dpts[l*DATAFIELD+3],dpts[l*DATAFIELD+4],dpts[l*DATAFIELD+1],dpts[l*DATAFIELD],
#else
			 p.qx,p.qy,p.qz,dpts[l*DATAFIELD+1],dpts[l*DATAFIELD],
#endif
			 results[l*(COUNTFIELD)+MU],results[l*(COUNTFIELD)+N],
			 results[l*(COUNTFIELD)+EP],results[l*(COUNTFIELD)+EN],results[l*(COUNTFIELD)+EM],results[l*(COUNTFIELD)+R],
			 results[l*(COUNTFIELD)+PRIMARYFIELD+S],results[l*(COUNTFIELD)+PRIMARYFIELD+G],
			 results[l*(COUNTFIELD)+PRIMARYFIELD+E]);
		}

	      fflush(stdout);

	      if (results)
		{
		  free(results);
		  results=NULL;
		}

#ifndef LUMP
	    }
	}
    }
#endif

  
  /* clean up */

  fprintf(stderr,"Cleaning up and exiting.\n");
  

  if (dpts)
    free(dpts);
  stop_nint();
  stop_ninterp();
  return EXIT_SUCCESS;
  
};












void * masterRHS( void * par, void * in, void * out )
{
  struct Parameters * p = par;
  double * retval = out;
  double * input = in;
  double smag;
  char * ecode = NULL;
  double lowerdomain;
#ifdef DEBUG
  fprintf(p->fout,"Made it to master. Finding mu...\n");
  fflush(p->fout);
#endif
  p->nc = input[1];
  p->J = input[0];
#ifdef LUMP
  p->qx = input[2];
  p->qy = input[3];
  p->qz = input[4];
#endif

  lowerdomain = -6. - 1.5*(p->J)*(p->magM);
  if ( lowerdomain > -10. )
    {
      lowerdomain = -10.;
    }
  retval[MU] = FINDMU(p,p->nc,lowerdomain,-1.*lowerdomain);
  p->mu = retval[MU];
  retval[N]  = quad( occupation_wrapper,p,lowerdomain,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
#ifdef DEBUG
      fprintf(p->fout,"Error in quad (occupation) : %s\n",ecode);
      fflush(p->fout);
#endif
      free(ecode);
      ecode=NULL;
    }
  p->nc = retval[N];
  /* paramagnetic calculations should be considered as broken. AFMISO code does not use occupation per-spin-channel, while
   * the functions in bare.c are for single electrons, [0,1]. For per-spin-channel code, the [0,1] range works fine. For
   * normal occupations, [0,2] and [0,4], I must put in place accounting for this, which I have not. So, in passing p->nc
   * on [0,4], paramagnetic_energy will be solving for nc on [0,1]. */ 
  //retval[EP] = paramagnetic_energy(p,LB,UB,&(retval[EN]),&(retval[EM]));
  retval[EP]=0.;  /* to be reinstated */ 
  retval[EN]=0.;    /* when paramagnetic */ 
  retval[EM]=0.;       /* energy is fixed */
  retval[R] = quad( iso_SCE_wrapper,p,lowerdomain,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
#ifdef DEBUG
      fprintf(p->fout,"Error in quad (iso_SCE) : %s\n",ecode);
      fflush(p->fout);
#endif
      free(ecode);
      ecode=NULL;
    }
  
  if ( retval[R] < 0. )
    {
      smag = p->J * (-1.*p->magM) * retval[R];
      p->magS = smag;
      retval[PRIMARYFIELD+S] = smag;
      retval[PRIMARYFIELD+G] = retval[MU]+1.5*p->J*smag; 
      p->g = retval[PRIMARYFIELD+G];
      retval[PRIMARYFIELD+E] = quad( eigen_occ_wrapper,p,lowerdomain,p->mu,0,REL_ERROR,&ecode);
      if (ecode)
	{
#ifdef DEBUG
      fprintf(p->fout,"Error in quad (eigen_occ): %s\n",ecode);
      fflush(p->fout);
#endif
	  free(ecode);
	  ecode=NULL;
	}
    }
#ifdef DEBUG
      fprintf(p->fout,"Finished in masterRHS.\n");
      fflush(p->fout);
#endif
  return (void *) retval;
}








void egress(int ac, ... )
{
  va_list alist;
  int i;
  int retval;
  void * container;
  if ( ac==0 )
    return;
  va_start(alist,ac);
  /* first list-argument must be the return value from dispatch */
  retval = va_arg(alist,int);
  /* all subsequent arguments are pointers to memory which must be freed */
  for (i=1;i<ac;i++)
    {
      container = (void *)va_arg(alist,void *);
      if ( container )
	{
	  free(container);
	  container=NULL;
	}
    }
  if ( retval == 0 )
    return;
  else
    {
      //finish_files();
      stop_nint();
      stop_ninterp();
      exit(retval);
    }
}





