/* Author: Gabriele Inghirami g.inghirami@gsi.de (2019-2022) - License: GPLv3 */

#ifndef DEF_MAIN
#include "definitions.h"
#endif

/** @file tools.c
* 
*   @brief this file contains the utilities, for example to initialize an array or to check if the files in the list given at the program invocation all exist
*/


extern const int start_index;

void check_input_files(char **infiles, int narg)
{
  int i;
  FILE *fp;
  for (i=start_index;i<narg;i++)
  {
    fp=fopen(infiles[i],"r");
    if(fp == NULL)
    {
      printf("File %s does not exist. I quit.\n",infiles[i]);
      exit(1);
    }
    else
    {
      fclose(fp);
    }
             
  }
}


void get_timesteps(char *timefile, int *ntimesteps, double **timesteps_array)
{
double *time_intervals_temporary_array;
const int size_of_time_intervals_temporary_array=100000;
int i;
FILE *ftint;

//we create a temporary array wich almost certainly is large enough to contain all the time intervals to be read
time_intervals_temporary_array=malloc(sizeof(double)*size_of_time_intervals_temporary_array);
ftint=fopen(timefile,"r");
printf("Opening timefile %s\n",timefile);
if( ftint== NULL)
{ 
  printf("Sorry, but I am unable to open the file with time intervals: %s\n",timefile);
  exit(1);
}
else
{
 *ntimesteps=0;
 while(fscanf(ftint,"%lf",&time_intervals_temporary_array[*ntimesteps]) != EOF) (*ntimesteps)++;
}
fclose(ftint);
*timesteps_array=malloc(sizeof(double)*(*ntimesteps));
for(i=0;i<*ntimesteps;i++) 
{
  (*timesteps_array)[i]=time_intervals_temporary_array[i];
}
//we delete the temporary array
free(time_intervals_temporary_array);
}


//it returns the index in the array of times time_int_array corresponding to the given test_time
int check_test_time(double test_time,double *time_int_array,int ntimesteps)
{
	int i;
	for(i=0;i<ntimesteps;i++)
	{
		if(time_int_array[i]==test_time) return i;
	}
	return -1;
}

