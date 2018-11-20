/*
  NAME:
     read_till_sep
  PURPOSE:
     reads a data value from the datafile
  CALLING SEQUENCE:
     read_till_sep(char curr_value[],FILE *file,char sep)
  INPUT:
     *file  - pointer to the datafile
     sep    - data separator in the datafile
  OUTPUT:
     curr_value - value that was read
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdbool.h>
#include <stdio.h>
#include <proj_gauss_main.h>

bool read_till_sep(char curr_value[],FILE *file,char sep){
  int vv=0;
  bool not_found_sep = true;
  char curr_char;
  while (not_found_sep){
    curr_char = (char) getc(file);
    if (curr_char == sep) break;
    if (curr_char == '\n') return true;
    curr_value[vv++] = curr_char;
  }

  return false;
}
