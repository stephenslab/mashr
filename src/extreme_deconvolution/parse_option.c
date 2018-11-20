/*
  NAME:
     parse_option
  PURPOSE:
     parse an option line in the initial conditions file
  CALLING SEQUENCE:
     parse_option(char line[])
  INPUT:
     line   - the line to parse
  OUTPUT:
     sets the relevant option
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdbool.h>
#include <string.h>
#include <proj_gauss_main.h>

bool parse_option(char line[]){
  //Define the options
  char K_opt[]="K";
  char fixamp_opt[]="fixamp";
  char fixmean_opt[]="fixmean";
  char fixcovar_opt[]="fixcovar";
  char maxiter_opt[]="maxiter";
  char tol_opt[]="tol";
  char likeonly_opt[]="likeonly";
  char w_opt[]="w";
  char splitnmerge_opt[]="splitnmerge";

  
  //Split the option in the name and value
  int ii=0,valuestart;
  char option[12]="", value[1000]="";
  while (line[ii] != '='){
    option[ii]= line[ii];
    ++ii;
  }
  valuestart=++ii;
  while (line[ii] != '\n'){
    value[ii-valuestart]= line[ii];
    ++ii;
  }
  
  //Set parameters
  if (strcmp(option,K_opt) == 0){
    K=atoi(value);
    //Based on this information, allocated some matrices/vectors
    fixampP = (bool *) malloc(K * sizeof (bool));
    fixmeanP = (bool *) malloc(K * sizeof (bool));
    fixcovarP = (bool *) malloc(K * sizeof (bool));
    gaussians = (struct gaussian *) malloc(K * sizeof (struct gaussian) );
    if (fixampP == NULL || fixmeanP == NULL || fixcovarP == NULL || gaussians == NULL){
      printf("Allocation of arrays failed, not enough free memory?\n");
      return -1;
    }
  }
  else if (strcmp(option,maxiter_opt) == 0) maxiter=atoll(value);
  else if (strcmp(option,tol_opt) == 0) tol=atof(value);
  else if (strcmp(option,w_opt) == 0) w=atof(value);
  else if (strcmp(option,likeonly_opt) == 0) likeonly=atoi(value);
  else if (strcmp(option,splitnmerge_opt) == 0) splitnmerge=atoi(value);
  else if (strcmp(option,fixamp_opt) == 0){
    for (ii=0; ii != K; ++ii){
      if (value[ii] == '0') *fixampP=false;
      else *fixampP=true;
      ++fixampP;
    }
    fixampP -= K;
  }
  else if (strcmp(option,fixmean_opt) == 0){
    for (ii=0; ii != K; ++ii){
      if (value[ii] == '0') *fixmeanP=false;
      else *fixmeanP=true;
      ++fixmeanP;
    }
    fixmeanP -= K;
  }
  else if (strcmp(option,fixcovar_opt) == 0){
    for (ii=0; ii != K; ++ii){
      if (value[ii] == '0') *fixcovarP=false;
      else *fixcovarP=true;
      ++fixcovarP;
    }
    fixcovarP -= K;
  }
  else return false;
  
  return true;
}
