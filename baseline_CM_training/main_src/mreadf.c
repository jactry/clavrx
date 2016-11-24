/*$Id$*/
/*****************************************************************************
 Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3

 NAME: mreadf.c (src)

 PURPOSE:
  This subroutine, which is callable from Fortran, reads in a
  user specified number of integer data records of a user specifed
  size given a byte offset from a binary file.

 DESCRIPTION:
  Input:
        filename: the name of the file to read from
	start: the byte offset value
	nbytes_per_rec: the size of each integer record
	nrec: the number of integer records to read
	buffer: integer pointer or array allocated to contain nrec
	        elements by the calling routine.
		
  Output:
        buffer: contains the integer data read in from file

 AUTHORS:
 Michael J. Pavolonis
 National Oceanic and Atmospheric Administration

 COPYRIGHT
 GEOCAT

 Copyright (C) 2006  Michael J. Pavolonis
 National Oceanic and Atmospheric Administration

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.
******************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

#define NFILE_MAX 200

typedef struct {
  FILE *fptr;
} VFILE_UNIT_TABLE;

VFILE_UNIT_TABLE file_unit_table[NFILE_MAX];

void mread_open_(const char *filename, int32_t *funit)
{
  if ((file_unit_table[(int)*funit].fptr = fopen(filename,"rb")) == NULL) {
      fprintf(stderr,"Could not open file1: %s\n", filename);
      exit(EXIT_FAILURE);
  }
}

void mread_close_(int32_t *funit)
{
  fclose(file_unit_table[(int)*funit].fptr);
}

void mreadf_int_o_(int32_t *funit,
                   int32_t *start,
	           int32_t *nbytes_per_rec,
	           int32_t *nrec,
	           int32_t *buffer)
{
  int32_t nrec_read;
  
  if (fseek(file_unit_table[(int)*funit].fptr, (int)*start, SEEK_SET) != 0) {
    fprintf(stderr,"Cannot point buffer to offset position in file unit %d\n",*funit);
    exit(EXIT_FAILURE);
  }
  
  nrec_read = fread((void *)buffer, (int)*nbytes_per_rec, (int)*nrec, file_unit_table[(int)*funit].fptr);
  if (nrec_read < (int)*nrec) {
    fprintf(stderr,"WARNING: The number of records read (%d) from unit %d is less than the number "
    "requested (%d)\n",nrec_read,*funit,*nrec);
  }
  
}

void mreadf_int_(const char *filename,
                 int32_t *start,
	         int32_t *nbytes_per_rec,
	         int32_t *nrec,
	         int32_t *nrec_read,
	         int32_t *buffer)
{
  FILE *fptr;
  
  if ((fptr = fopen(filename,"rb")) == NULL) {
      fprintf(stderr,"Could not open file2: %s\n", filename);
      exit(EXIT_FAILURE);
  }
  
  if (fseek(fptr, (int)*start, SEEK_SET) != 0) {
    fprintf(stderr,"Cannot point buffer to offset position\n");
    exit(EXIT_FAILURE);
  }
  
 *nrec_read = fread(buffer, (int)*nbytes_per_rec, (int)*nrec, fptr);
  if (*nrec_read < (int)*nrec) {
    fprintf(stderr,"WARNING: The number of records read (%d) from %s is less than the number "
    "requested (%d)\n",*nrec_read,filename,*nrec);
  }
  
  fclose(fptr);
  
}
