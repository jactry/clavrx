/*$Id: move_bytes.c,v 1.4.2.2 2014/01/26 04:49:37 heidinger Exp $*/
/*****************************************************************************
 Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 5.3

 NAME: move_bytes.c (src)

 PURPOSE:

 DESCRIPTION:

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
#include<string.h>

void move_bytes_(int32_t *nbytes, void *inbuffer, void *outbuffer, int32_t *offset)
{
  unsigned char  *outbuf = outbuffer;
  memcpy(&outbuf[(int)*offset], inbuffer, (int)*nbytes);
}

void swap_bytes4_(void *buffer, int32_t *num)
{
  char *cbuf = (char *)buffer;
  int i, n;

  n = (int)*num;
  for (i = 0; i < n; i++)
  {
    char b;

    b = cbuf[0];
    cbuf[0] = cbuf[3];
    cbuf[3] = b;
    b = cbuf[1];
    cbuf[1] = cbuf[2];
    cbuf[2] = b;
    cbuf += 4;
  }
}

void swap_bytes2_(void *buffer, int32_t *num)
{
  char *cbuf = (char *)buffer;
  int i, n;

  n = (int)*num;
  for (i = 0; i < n; i++)
  {
    char b;

    b = cbuf[0];
    cbuf[0] = cbuf[1];
    cbuf[1] = b;
    cbuf += 2;
  }
}
