/**************************************************************************************************************
 *  fio.h - I/O utility library.
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
 ******************************************************************************************************************
 *
 *  Utility library.
 *  Handles multiple-file inputs:
 *  Local variables hold the file list so that each call to get_data operates on a single file.
 *  On get_data call returns all the data from one file.
 *  A subsequent get_data calls reads from the -next- file in the list.
 *  When the list has been exhausted, get_data returns 0 and a NULL pointer.
 *  By contrast, get_mesh operates with a single mesh file and returns only this information.
 *  To use, prep_file_list must be called before get_data,
 *    and finish_files must be called when the library calls are no longer needed.
 *  This library cannot handle the case where:
 *    a file is missing AFTER the file pointer has been placed in the file list
 *  In retrospect, this library would have been more easily and more cleanly implemented using an altered
 *    list.h from the linux headers.
 *  1 oct 2016
 **************************************************************************************************************/


#ifndef F_IO


#define F_IO

#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>


/* for multiple data files */
int prep_file_list( const char * d );


/* opens file "<d>/<f>" where d is the absolute directory path and f is the filename.
   Works independently of prep_file_list and finish_files */
int get_mesh( const char * d, const char * f, double ** data );


/* utility for mapping mesh data retrieved via get_mesh into two arrays (x,y).
   Need a more general version of this file to handle n-dimensional mesh data. */ 
int mesh_copy(double ** x, double ** y, double * dat, int n);


/* from the data files prepared in prep_file_list, retrieves the full data of a single file once per call.
   Automatically iterates through files, returns 0 and NULL when the list is exhausted.
   On a successful call, returns the number of records read. */
int get_data( double ** data );


/* closes data files associated with prep_file_list */
int finish_files( );


/* utility for opening an arbitrary list of files. Not related to mesh or data fctns.
   "n" must be the number of file names and number of FILE * spots. */
int open_files(FILE ** files, const char * const * names, int n);


/* utility for resizing an array. v is the address of the array pointer. sz is the size of each datum.
   ncur is the current length of v. return value is the new number of array elements. */
int resize_g( void ** v, const int sz, const int ncur );


#endif /* F_IO */
