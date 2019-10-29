/*******************************************************************************************************************************
 *  
 *  General zero-T K-PM state solver (zkpms)
 *  Altered 2018 june 13 as per PSR's argument for sign of linearization constant, etc.
 *
 *  2018 jun 6
 *    o  adapted code to "lump" data...for owl's nest. Build with -DLUMP to enable.
 *  2018 mar 1
 *    o  removed sce-return-val from output; replaced with free energy constant (C) and eigen-occup val (EIG)
 *    o  applied changes to integration method as well (not done previously)
 *    o  output width is changed to 1.6e from 1.8e
 *
 *******************************************************************************************************************************/


#include <stdio.h>
#include <stdarg.h>
#include "nint.h"
#include "fio.h"
#include "sm.h"
#include "kzerocore.h"

/* unistd.h is only used in this file if debugging is needed */
#ifdef DEBUG
#  include <sys/unistd.h>
#endif


#define LOWERJV           ( 0.002 )  /* J>0, V>0 */
#define INCJV             ( .001 )
#define COUNTJV           ( 1 )
#define LOWERG            ( 0. )  /* Increase the bounds of g to achieve higher & lower nc's */
#define COUNTG            ( 1 )
#define INCG              ( .2 )
#define COUNTFIELDS       ( 22 )
#ifdef LUMP
#  define DATAFIELD       ( 2 )
#else
#  define DATAFIELD       ( 1 )
#endif

const char dirname[8] = "./input\0";
const char meshname[17] = "dos_mesh_t_1.txt\0";


void * masterSCEV ( void * par, void * in, void * out );
void dummyhook(int i);
void egress(int, ... );



enum {           /* return value indices; L -> lower possible mu ; U -> upper possible mu if applicable */
  LJ,            /* J */
  LV,            /* V */
  LMU,           /* mu */
  LNF,           /* nf */
  LNC,           /* nc */
  LE0,           /* PM phase energy */
  LEN,           /* PM nc for Helmholtz calculation */
  LEM,           /* PM mu for Helmholtz calculation */
  LC,            /* Constants in the Free energy */
  LEIG,          /* Eigen-occup-c return value */
  LE,            /* Helmholtz Free energy at T=0 */
  UJ,
  UV,
  UMU,
  UNF,
  UNC,
  UE0,
  UEN,
  UEM,
  UC,
  UEIG,
  UE
};



int
main (void)
{

  /* helper vars */
  int i = 0,j=0,k;
  int count, mcount, jvcount; /* total numbers of input test data, mu values, and JV values, resp. */
  int ret;
  char * errorcode;           /* for GSL error-codes passed through "quad" function */

  double * dpts = NULL;
  double * results = NULL;
  double * meshx, * meshy;
  int meshn;
  int datacount;
  double jv ;
  double sce_ret;
  double E;
  double Enot;
  
  struct Parameters p;
  struct d_Function df; /* for dispatching, from sm library */

#ifdef DEBUG
  p.fout = stdout;
#endif


  /* fetch the d.o.s. interpolation mesh */

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

  
  /* prepare numerical integration */
  
  start_ninterp(meshn,meshx,meshy);
  meshx=NULL;
  meshy=NULL;
  start_nint(NINTINT);


  /* diagnostics on occupation */
  /*
  double inc=-6.;
  p.J=-3.58157144;
  p.V=.293167404;
  p.g=-0.1;
  while ( inc<6. )
    {
      p.mu=inc;
#ifdef KSUM
      //printf("%f %1.6e\n",inc,SUMMOR(occupation_bare_zero_k,&p,SUMPREC));
      printf("%f %1.6e\n",inc,k_trapezoidal(eigen_occup_k,&p,.1));
#else
      //printf("%f %1.6e %1.6e\n",inc,quad(rho_c,&p,LOWERINT,inc,0,REL_ERROR,&errorcode),quad(rho_m,&p,LOWERINT,inc,0,REL_ERROR,&errorcode));
      printf("%f %1.6e\n",inc,quad(occupation_bare_zero,NULL,LB,inc,0,REL_ERROR,&errorcode));
#endif
      inc += .05;
    }
#ifndef KSUM
  stop_nint();
  stop_ninterp();
#endif
  return EXIT_SUCCESS;*/
  

  /* informational printout */

#ifdef USECASEJ32
  printf("#J=3/2 ");
  fprintf(stderr,"J=3/2 ");
#  ifdef USEMAGCASE32
  printf(" (Mag 3/2 is defined but useless here)\n");
  fprintf(stderr," (Mag 3/2 is defined but useless here)\n");
#  else /*notmagcase32*/
  printf("\n");
  fprintf(stderr,"\n");
#  endif
#else
  printf("#S=1/2\n");
  fprintf(stderr,"S=1/2\n");
#endif
  printf("#Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  fprintf(stderr,"Precision: %1.2e; k-sum precision: %1.2e; relative error: %1.2e; nint intervals: %d\n",PREC,SUMPREC,REL_ERROR,NINTINT);
  printf("#NIMP is %f\n",NIMP);
  fprintf(stderr,"NIMP is %f\n",NIMP);
  printf("#%d(JV) x %d(lambda) = %d(records)\n",COUNTJV,COUNTG,COUNTJV*COUNTG);
  fprintf(stderr,"%d(JV) x %d(lambda) = %d(records)\n",COUNTJV,COUNTG,COUNTJV*COUNTG);
  

  /* print banner */

  printf("##g    JV     mu           nf           nc           J             V            E0            EN           EM           C             EIG          E\n");

  
  /* prepare parameter space */

  fprintf(stderr,"Creating parameter space...\n");
#ifdef LUMP
  datacount = COUNTG*COUNTJV;
#else
  datacount = COUNTJV;
#endif
  dpts = malloc(sizeof(double)*DATAFIELD*datacount);
  if ( dpts==NULL )
    {
      printf("Failed malloc for dpts\n");
      egress(1,-1);
    }
#ifdef LUMP
  k = 0;
  for (i=0;i<COUNTJV;i++)
    {
      for (j=0;j<COUNTG;j++)
	{
	  dpts[k++] = (LOWERJV+((double)i)*INCJV) ;
	  dpts[k++] = (LOWERG+((double)j)*INCG) ;
	}
    }
#else
  for (i=0;i<COUNTJV;i++)
    {
      dpts[i] = (LOWERJV+((double)i)*INCJV) ;
    }
#endif


  /* prepare sm (dispatching) data structures */

  df.f = masterSCEV ;
  df.param = &p;


  /* dispatch loop */

  fprintf(stderr,"Looping\n");
  fflush(stdout);
  fflush(stderr);

#ifndef LUMP
  for (j=0;j<COUNTG;j++)
    {

      p.g = LOWERG + ((double)j)*INCG;
#endif
      
      dispatchv2(df,datacount,dpts,sizeof(double)*DATAFIELD,&results,COUNTFIELDS*sizeof(double),egress);

      
      for (i=0;i<datacount;i++)
	{
	  /*      g     jv    lmu   lnf   lnc   lJ    LV    lE0   lEN   lEMu  LC    LEIG  lE                                      */	  
	  printf("%1.2f %1.4f %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e\n",
#ifdef LUMP
		 dpts[i*DATAFIELD+1],dpts[i*DATAFIELD],
#else
		 p.g,dpts[i],
#endif
		 results[i*COUNTFIELDS+LMU],results[i*COUNTFIELDS+LNF],results[i*COUNTFIELDS+LNC],
		 results[i*COUNTFIELDS+LJ],results[i*COUNTFIELDS+LV],results[i*COUNTFIELDS+LE0],results[i*COUNTFIELDS+LEN],
		 results[i*COUNTFIELDS+LEM],results[i*COUNTFIELDS+LC],results[i*COUNTFIELDS+LEIG],results[i*COUNTFIELDS+LE]);
	  if (results[LMU]!=results[UMU])
	    {
		 printf("%1.2f %1.4f %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6e %1.6 %1.6e\n",
#ifdef LUMP
		        dpts[i*DATAFIELD+1],dpts[i*DATAFIELD],
#else
		        p.g,dpts[i],
#endif
		        results[i*COUNTFIELDS+UMU],results[i*COUNTFIELDS+UNF],results[i*COUNTFIELDS+UNC],
		        results[i*COUNTFIELDS+UJ],results[i*COUNTFIELDS+UV],results[i*COUNTFIELDS+UE0],results[i*COUNTFIELDS+UEN],
		        results[i*COUNTFIELDS+UEM],results[i*COUNTFIELDS+UC],results[i*COUNTFIELDS+UEIG],results[i*COUNTFIELDS+UE]);
	    }
	}

      fflush(stdout);

      if (results)
	{
	  free(results);
	  results=NULL;
	}

#ifndef LUMP
    }
#endif

  
  /* clean up */

  if (dpts)
    {
      free(dpts);
    }
  if (results)
    {
      free(results);
    }
  stop_nint();
  stop_ninterp();
  
  return EXIT_SUCCESS;
}




#ifdef KSUM


void * masterSCEV ( void * par, void * dat, void * rdat )
{
  struct Parameters * p = par;
  double * input = dat;
  double * output = rdat;
  double uppermu,lowermu,sce_ret;
  p->J = input[0];
#ifdef LUMP
  p->g = input[1];
#endif
  p->V = 1.;
  if (msplit_k(p,LB,UB,&lowermu,&uppermu)<0)
    {
      /* error in GSL */
      output[LMU] = 2*LOWER;
      output[LJ] = 0.;
      output[LV] = 0.;
      output[LNF] = 0.;
      output[LNC] = 0.;
      output[LE0] = 0.;
      output[LEN] = 0.;
      output[LEM] = 0.;
      output[LC] = 0.;
      output[LEIG] = 0.;
      output[LE] = 0.;      
      output[UMU] = 2*LOWER;
      output[UJ] = 0.;
      output[UV] = 0.;
      output[UNF] = 0.;
      output[UNC] = 0.;
      output[UE0] = 0.;
      output[UEN] = 0.;
      output[UEM] = 0.;
      output[UC] = 0.;
      output[UEIG] = 0.;
      output[UE] = 0.;
      return (void *)output;
    }
  if ( lowermu>UPPER )
    {
      output[LMU] = lowermu;
      output[LJ] = 0.;
      output[LV] = 0.;
      output[LNF] = 0.;
      output[LNC] = 0.;
      output[LE0] = 0.;
      output[LEN] = 0.;
      output[LEM] = 0.;
      output[LC] = 0.;
      output[LEIG] = 0.;
      output[LE] = 0.;      
      output[UMU] = lowermu;
      output[UJ] = 0.;
      output[UV] = 0.;
      output[UNF] = 0.;
      output[UNC] = 0.;
      output[UE0] = 0.;
      output[UEN] = 0.;
      output[UEM] = 0.;
      output[UC] = 0.;
      output[UEIG] = 0.;
      output[UE] = 0.;
      return (void *)output;
    }
  p->mu = lowermu;
  output[LMU]=lowermu;
  output[LNC]=SUMMOR(rho_c_k,p,SUMPREC);
  p->nc = output[LNC];
  output[LNF]=SUMMOR(rho_m_k,p,SUMPREC);
  sce_ret=SUMMOR(sce_k,p,SUMPREC);
  p->J = 1./sce_ret;
  p->V = sce_ret*(*input);
  output[LJ]=p->J;
  output[LV]=p->V;
  output[LE0] = paramagnetic_energy(p,LB,UB,&(output[LEN]),&(output[LEM]));
#  ifdef USECASEJ32
  output[LC] = 3.75*(*input)*(p->V);
  output[LEIG] = SUMMOR(eigen_occup_k,p,SUMPREC); 
  output[LE]= output[LEIG] - .25*(p->g) + output[LC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  else /*S=1/2*/
  output[LC] = .75*(*input)*(p->V);
  output[LEIG] = SUMMOR(eigen_occup_k,p,SUMPREC);
  output[LE]= output[LEIG] - .5*(p->g) + output[LC] ;  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  endif /*USECASEJ32*/
  if ( uppermu==lowermu )
    {
      output[UMU]=lowermu;
      output[UJ]=0.;
      output[UV]=0.;
      output[UNF]=0.;
      output[UNC]=0.;
      output[UE0]=0.;
      output[UEN]=0.;
      output[UEM]=0.;
      output[UC]=0.;
      output[UEIG]=0.;
      output[UE]=0.;
      return (void *)output;
    }
  p->mu = uppermu;
  output[UMU] = uppermu;
  output[UNC]=SUMMOR(rho_c_k,p,SUMPREC);
  p->nc = output[UNC];
  output[UNF]=SUMMOR(rho_m_k,p,SUMPREC);
  sce_ret=SUMMOR(sce_k,p,SUMPREC);
  p->J = 1./sce_ret;
  p->V = sce_ret*(*input);
  output[UJ]=p->J;
  output[UV]=p->V;
  output[UE0] = paramagnetic_energy(p,LB,UB,&(output[UEN]),&(output[UEM]));
#  ifdef USECASEJ32
  output[UC] = 3.75*(*input)*(p->V);
  output[UEIG] = SUMMOR(eigen_occup_k,p,SUMPREC); 
  output[UE]= output[UEIG] - .25*(p->g) + output[UC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  else /*S=1/2*/
  output[UC] = .75*(*input)*(p->V);
  output[UEIG] = SUMMOR(eigen_occup_k,p,SUMPREC);
  output[UE]= output[UEIG] - .5*(p->g) + output[UC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  endif /*USECASEJ32*/
  return (void *)output;
}



#else /*not KSUM*/



void * masterSCEV ( void * par, void * dat, void * rdat )
{
  struct Parameters * p = par;
  double * input = dat;
  double * output = rdat;
  double uppermu,lowermu,sce_ret;
  char * ecode = NULL;
  p->J = input[0];
#ifdef LUMP
  p->g = input[1];
#endif
  p->V = 1.;
  if (msplit(p,LB,UB,&lowermu,&uppermu)<0)
    {
      /* error in GSL */
      output[LMU] = 2*LOWER;
      output[LJ] = 0.;
      output[LV] = 0.;
      output[LNF] = 0.;
      output[LNC] = 0.;
      output[LE0] = 0.;
      output[LEN] = 0.;
      output[LEM] = 0.;
      output[LC] = 0.;
      output[LEIG] = 0.;
      output[LE] = 0.;      
      output[UMU] = 2*LOWER;
      output[UJ] = 0.;
      output[UV] = 0.;
      output[UNF] = 0.;
      output[UNC] = 0.;
      output[UE0] = 0.;
      output[UEN] = 0.;
      output[UEM] = 0.;
      output[UC] = 0.;
      output[UEIG] = 0.;
      output[UE] = 0.;
      return (void *)output;
    }
  if ( lowermu>UPPER )
    {
      output[LMU] = lowermu;
      output[LJ] = 0.;
      output[LV] = 0.;
      output[LNF] = 0.;
      output[LNC] = 0.;
      output[LE0] = 0.;
      output[LEN] = 0.;
      output[LEM] = 0.;
      output[LC] = 0.;
      output[LEIG] = 0.;
      output[LE] = 0.;      
      output[UMU] = lowermu;
      output[UJ] = 0.;
      output[UV] = 0.;
      output[UNF] = 0.;
      output[UNC] = 0.;
      output[UE0] = 0.;
      output[UEN] = 0.;
      output[UEM] = 0.;
      output[UC] = 0.;
      output[UEIG] = 0.;
      output[UE] = 0.;
      return (void *)output;
    }
  p->mu = lowermu;
  output[LMU]=lowermu;
  output[LNC]=quad(rho_c,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  p->nc = output[LNC];
  output[LNF]=quad(rho_m,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  sce_ret=quad(sceV,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  p->J = 1./sce_ret;
  p->V = sce_ret*(*input);
  output[LJ]=p->J;
  output[LV]=p->V;
  output[LE0] = paramagnetic_energy(p,LB,UB,&(output[LEN]),&(output[LEM]));
#  ifdef USECASEJ32
  output[LC] = 3.75*(*input)*(p->V);
  output[LEIG] = 0.5 * quad(eigen_occup,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  output[LE]= output[LEIG] - .25*(p->g) + output[LC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  else /*S=1/2*/
  output[LC] = .75*(*input)*(p->V);
  output[LEIG] = 0.5 * quad(eigen_occup,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  output[LE]= output[LEIG] - .5*(p->g) + output[LC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  endif /*USECASEJ32*/
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  if ( uppermu==lowermu )
    {
      output[UMU]=lowermu;
      output[UJ]=0.;
      output[UV]=0.;
      output[UNF]=0.;
      output[UNC]=0.;
      output[UE0]=0.;
      output[UEN] = 0.;
      output[UEM] = 0.;
      output[UC] = 0.;
      output[UEIG] = 0.;
      output[UE]=0.;
      return (void *)output;
    }
  p->mu = uppermu;
  output[UMU] = uppermu;
  output[UNC]=quad(rho_c,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  p->nc = output[UNC];
  output[UNF]=quad(rho_m,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  sce_ret=quad(sceV,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  p->J = 1./sce_ret;
  p->V = sce_ret*(*input);
  output[UJ]=p->J;
  output[UV]=p->V;
  output[UE0] = paramagnetic_energy(p,LB,UB,&(output[UEN]),&(output[UEM]));
#  ifdef USECASEJ32
  output[UC] = 3.75*(*input)*(p->V);
  output[UEIG] = 0.5 * quad(eigen_occup,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  output[UE]= output[UEIG] - .25*(p->g) + output[UC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  else /*S=1/2*/
  output[UC] = .75*(*input)*(p->V);
  output[UEIG] = 0.5 * quad(eigen_occup,p,LOWERINT,p->mu,0,REL_ERROR,&ecode);
  output[UE]= output[UEIG] - .5*(p->g) + output[UC];  // 4 (15/4) JV^2 + 4 eigen*occup ;
#  endif /*USECASEJ32*/
  if (ecode)
    {
      free(ecode);
      ecode==NULL;
    }
  return (void *)output;
}


#endif /*KSUM*/







void dummyhook(int i)
{
  return;
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
	free(container);
    }
  if ( retval == 0 )
    return;
  else
    {
      finish_files();
      stop_nint();
      stop_ninterp();
      exit(retval);
    }
}


