/**************************************************************************************************************
 *  fio.c - I/O utility library.
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
 ******************************************************************************************************************/

#include "fio.h"


#define START_SZ 6
#define START_DATA_SZ 512
#define PATH_SZ 256
int f_list_pos;
int f_n_tot;
FILE ** f_list;


struct bundle {
  double a,b,c;
};


/* resize file arrays. Just use resize_g instead */
int resize( FILE *** farray, const int cursize)
{
  int newsize = 2*cursize + 1;
  int i;
  FILE ** retarray = malloc(sizeof(FILE *)*newsize);
  if ( retarray==NULL )
    {
      return -1;
    }
  for (i=0;i<cursize;i++)
    {
      retarray[i] = (*farray)[i];
    }
  for (i;i<newsize;i++)
    {
      retarray[i] = NULL;
    }
  free(*farray);
  *farray = retarray;
  return newsize;
}


/* general resize: v is the array of pointers to be resized, sz is the size of each element in the array, ncur is the current number of elements */
/* return value is the new number of array elements */
int resize_g( void ** v, const int sz, const int ncur )
{
  int newsize = 2*ncur + 1;
  int i;
  void * retarray = malloc(sz*newsize);
  if ( retarray==NULL )
    {
      return -1;
    }
  bcopy( *v, retarray, sz*ncur );
  bzero( (void *)(retarray+(ncur*sz)), sz*ncur+1 );
  free(*v);
  *v = retarray;
  return newsize;
}


int parse_name( const char * name, int n )
{
  /* regex "nn*.txt" */
  char * r;
  r = strstr(name,"nn");
  if ( r==name )
    {
      r = strstr(name,".txt");
      if (r == &(name[n-4]))
	{
	  return 1;
	}
    }
  return 0;
}


int prep_file_list( const char * d )
{
  struct dirent * file = NULL;
  int fcounter = 0;
  int fcursz = START_SZ;
  char fullpath[PATH_SZ];
  DIR * dir = NULL;
  
  f_list = NULL;
  f_list_pos = 0;
  dir = opendir(d);
  if (dir==NULL)
    {
      return -1;
    }
  f_list = malloc(sizeof(FILE *)*fcursz);
  if ( f_list==NULL )
    {
      closedir(dir);
      return -2;
    }
  bzero(f_list,sizeof(FILE *)*fcursz);
  bzero(fullpath,PATH_SZ*sizeof(char));
  while ( (file=readdir(dir), file)!=NULL )
    {
      if (parse_name(file->d_name,strlen(file->d_name))>0)
	{
	  sprintf(fullpath,"%s/%s",d,file->d_name);
	  f_list[fcounter]=fopen(fullpath,"r");
	  fcounter++;
	  if ( fcounter==fcursz )
	    {
	      fcursz = resize_g((void **)&f_list,sizeof(FILE *),fcursz);
	      if ( fcursz<1 )
		{
		  closedir(dir);
		  for (fcounter;fcounter>=0;fcounter--)
		    {
		      free(f_list[fcounter]);
		    }
		  free(f_list);
		  return -4;
		}
	    }
	}
      else
	continue;
    }
  f_n_tot = fcounter;
  closedir(dir);
  return fcounter;
}


int get_mesh( const char * d, const char * f, double ** data )
{
  char fullpath[PATH_SZ];
  FILE * fmesh;
  int dcursz = START_DATA_SZ*2;
  int i;
  
  bzero(fullpath,PATH_SZ*sizeof(char));
  sprintf(fullpath,"%s/%s",d,f);
  fmesh = fopen(fullpath,"r");
  if (fmesh==NULL)
    return -1;
  *data = malloc(dcursz*sizeof(double));
  if ( *data==NULL )
    {
      fclose(fmesh);
      return -2;
    }
  i=-1;
  while (feof(fmesh)!=EOF)
    {
      i++;
      if ( ((i+1)*2)>=dcursz )
	dcursz = resize_g( (void **)(data), sizeof(double), dcursz );
      if ( fscanf(fmesh,"%lf %lf\n",&((*data)[i*2]),&((*data)[i*2+1])) < 1 )
	{
	  break;
	}
    }
  if ( ferror(fmesh) )
    {
      fclose(fmesh);
      return -3;
    }
  fclose(fmesh);
  return i;
}


int
mesh_copy(double ** x, double ** y, double * dat, int n)
{
  int i;
  *x = malloc( n*sizeof(double) );
  if ( *x==NULL )
    {
      return -1;
    }
  *y = malloc( n*sizeof(double) );
  if ( *y==NULL )
    {
      free(*x);
      return -1;
    }
  for (i=0;i<n;i++)
    {
      (*x)[i] = dat[2*i];
      (*y)[i] = dat[2*i+1];
    }
  return 0;
}


int get_data( double ** data )
{
  int i;
  int dcounter = 0;
  int dcursz = START_DATA_SZ*3;
  if ( data == NULL )
    {
      return -3;
    }
  if ( f_list_pos >= f_n_tot )
    {
      *data = NULL;
      return 0;
    }
  *data = malloc(dcursz*sizeof(double));
  if ( *data==NULL )
    {
      return -1;
    }
  i=-1;
  while (feof(f_list[f_list_pos])!=EOF)
    {
      i++;
      if ( ((i+1)*3)>=dcursz )
	dcursz = resize_g( (void **)(data), sizeof(double), dcursz );
      if ( fscanf(f_list[f_list_pos],"%lf %lf %lf\n",&((*data)[i*3]),&((*data)[i*3+1]),&((*data)[i*3+2])) < 1 )
	{
	  break;
	}
    }
  if ( ferror(f_list[f_list_pos++]) )
    return -2;
  return i;
}


int finish_files( )
{
  int i;
  if ( f_list )
    {
      for (i=0;i<f_n_tot;i++)
	if ( f_list[i] )
	  fclose(f_list[i]);
      free(f_list);
    }
  return 0;
} 


int
open_files(FILE ** files, const char * const * names, int n)
{
  int j = 0;
  int k = 0;
  do
    {
      files[j]=fopen(names[j],"w");
      if (files[j]==NULL)
	{
	  while (k<3)
	    {
	      if (files[k]!=NULL)
		fclose(files[k]);
	      else
		break;
	    }
	  return -1;
	}
      j++;
    } while (j<n);
  return 0;
}
