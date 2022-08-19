#ifndef DEF_MAIN
#include "definitions.h"
#endif
/**
* @file main.c
*
* @brief this file contains the main function, which controls the highest level of the program execution. It also contains the help function. 
*/

#if USE_CENTERED_GRID==0
double const xmin=-NX_DEF*DX_DEF/2.;
double const ymin=-NY_DEF*DY_DEF/2.;
double const zmin=-NZ_DEF*DZ_DEF/2.;
#else
double const xmin=X_START;
double const ymin=Y_START;
double const zmin=Z_START;
#endif

double const dx=DX_DEF;
double const dy=DY_DEF;
double const dz=DZ_DEF;
int const nx = NX_DEF;
int const ny = NY_DEF;
int const nz = NZ_DEF;
int const np = NP+1;

const char Tplabel[]="_Tmunu_";
const char infolabel[]="_Tmunu_info.dat";
const char dens_label[]="_densities_";
const char info_dens_label[]="_densities_info.dat";

int nt;
size_t pdata_totmem=0;
size_t big_arrays_allocated_mem=0;
size_t max_memory_allocatable_data;

const int T10=T01;
const int T20=T02;
const int T30=T03;
const int T21=T12;
const int T31=T13;
const int T32=T23;

short const b_selection=B_SELECTION;
float_t const bmin = BMIN;
float_t const bmax = BMAX;

const double cell_volume=DX_DEF*DY_DEF*DZ_DEF;

int const max_args = 2009; //maximum number of arguments 

const int index_of_output_file=3;//index of the output file in argv

const int start_index=4;//index of the argument reporting the first input file. 

const int shift_resonance_index = 10000; //index offset to distinguish resonances and stable particles

double *time_int_array; //array with the selected times

/** @brief main function
* It parses the command line arguments, it allocates the most important data arrays (Tp, Jp, Jb, Jc and Js) and it calls the compute or avg functions.
* @param[in] <comp|avg> comp : it computes the energy momentum tensor for all particles, the baryon and the individual particle four currents;  avg : it averages previously computed datafiles
* @param[in] <0|1> 0: It does not print separate baryon density current files 1: it does
* @param[in] <outpufile_prefix> The root of the names of the ouput files (the others part of their names are built automatically)
* @param[in] <UrQMD_event_file_1|Tmunu_file1> The first UrQMD event file (.f14) if the program is launched with the comp option or the first T_munu output file to be averaged if the program is launched with the avg option
* @param[in] <UrQMD_event_file_2|Tmunu_file2> The second UrQMD event file (.f14) if the program is launched with the comp option or the second T_munu output file to be averaged if the program is launched with the avg option. The program accepts up to roughly 1000 events files, but the OS might impose a tighter limit.
* @param[in] <file_with_time_intervals|Tmunu_fileN> Either the file with the chosen time intervals to compute ( program launched with comp ) or the final T_munu file to be averaged ( program launched with avg )
* \callgraph
**/

int main(int argc, char *argv[])
{


double *Tp; //< The array containing the T_munu tensor for all particle species. It is a linear array of dimensions: nt*nx*ny*nz*np*10, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively, np is the number of particles and 10 corresponds to the ten energy momentum tensor components for each particle. 
double *Jp; //< The array containing the four current of all particle species. It is a linear array of dimensions: nt*nx*ny*nz*np*4, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively, np is the number of particles and 4 corresponds to the four covariant components for each particle.
double *Jb; //< The array containing the net baryon four current. It is a linear array of dimensions: nt*nx*ny*nz*4, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively and 4 corresponds to the four covariant components.
double *Jc; //< The array containing the net electric four current. It is a linear array of dimensions: nt*nx*ny*nz*4, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively and 4 corresponds to the four covariant components.
double *Js; //< The array containing the net strangeness four current. It is a linear array of dimensions: nt*nx*ny*nz*4, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively and 4 corresponds to the four covariant components.
double *Jt; //< The array containing the total baryon four current. It is a linear array of dimensions: nt*nx*ny*nz*4, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively and 4 corresponds to the four covariant components.
long int *Pnum; //< The array containing the total number of individual particle species. It is a linear array of dimensions: nt*nx*ny*nz*np, where nt is the number of timesteps, nx, ny and nz are the cells along x, y and z respectively, np is the number of particles. Each entry tells the total number of particles of the species p.

long int nevents=0;//< The number of events
char * outputfile;
int print_dens;//< The flag to decide whether to print or not the four current components into additional separate files ( the four currents are always saved with the energy momentum tensors, but these additional files make easier to extract their values in text format )


//in the worst case, the number of arguments must be at least 5
if(argc<5) 
{
  help();
  exit(1);
}

if(argc > max_args)
{
  printf("Sorry, but opening too many files can cause problems, so there is a limit (%d) on the number of accepted arguments...\nIf you wish, you can modify the constant int max_args, recompile the program and run it again.\n",max_args);
  return 1;
}

if(strlen(argv[3])>1000) //this check might be looser or removed, but the length should be already sufficient for most purposes
{
	printf("The name of the output prefix cannot exceed 1000 characters.\n");
	exit(1);
}
else
{
	outputfile=argv[3];
}


max_memory_allocatable_data=sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGE_SIZE)*MAX_MEMORY_ALLOC;

//formal check of the parameters
if((strncmp(argv[1],"comp",4) == 0) || (strncmp(argv[1],"avg",3) == 0))
{
  check_input_files(argv,argc);//we check that the input files really exist, the final numer tells from which argument index to start

  if(strncmp(argv[1],"comp",4) == 0) {
    get_timesteps(argv[argc-1],&nt,&time_int_array);
  }
  else {
     time_int_array=(double *)calloc(1,sizeof(double));
     nt=1;
  }

  Tp=(double *)calloc(nt*nx*ny*nz*np*10,sizeof(double));
  if(Tp==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Tp array inside main. I am forced to quit.\n");
	  exit(4);
  }
  Jp=(double *)calloc(nt*nx*ny*nz*np*4,sizeof(double));
  if(Jp==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Jp array inside main. I am forced to quit.\n");
	  exit(4);
  }
  Jb=(double *)calloc(nt*nx*ny*nz*4,sizeof(double));
  if(Jb==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Jb array inside main. I am forced to quit.\n");
	  exit(4);
  }
  Jc=(double *)calloc(nt*nx*ny*nz*4,sizeof(double));
  if(Jc==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Jc array inside main. I am forced to quit.\n");
	  exit(4);
  }
  Js=(double *)calloc(nt*nx*ny*nz*4,sizeof(double));
  if(Js==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Js array inside main. I am forced to quit.\n");
	  exit(4);
  }
  Jt=(double *)calloc(nt*nx*ny*nz*4,sizeof(double));
  if(Jt==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Jt array inside main. I am forced to quit.\n");
	  exit(4);
  }
  Pnum=(long int *)calloc(nt*nx*ny*nz*np,sizeof(long int));
  if(Pnum==NULL)
  {
	  printf("Sorry, but it is not possible to allocate the Pnum array inside main. I am forced to quit.\n");
	  exit(4);
  }


  big_arrays_allocated_mem=((size_t)(14*sizeof(double)*nt*np)+(size_t)(4*sizeof(double)*nt))*((size_t)(nx*ny*nz));
  if(big_arrays_allocated_mem>max_memory_allocatable_data)
  {
	  printf("The memory required for the allocation of the main arrays is %zu.\n",big_arrays_allocated_mem);
	  printf("The total ram of this system is estimated to be %zu.\n",sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGE_SIZE));
	  printf("The ram memory usable by the program according to the parameter MAX_MEMORY_ALLOC in particle_definitions.h (which sets the fraction of the available memory) is %zu\n.",max_memory_allocatable_data);
	  printf("Please, check your settings again. I quit.\n");
	  exit(4);
  } 
}
else
{
  printf("Unknown option %s\n",argv[1]);
  help();
  exit(1);
}


if((strncmp(argv[1],"comp",4) == 0))
{  
  #ifdef SMASH
  prepare_smash_hadron_array(); 
  #endif
  compute(argv,argc,Tp,Jp,Jb,Jc,Js,Jt,Pnum,&nevents);  
}
else if(strncmp(argv[1],"avg",3) == 0)
{
  if(big_arrays_allocated_mem*2>max_memory_allocatable_data)
  {
	  printf("The memory required for the allocation of the main arrays (in the average procedure we need them twice) is %zu.\n",big_arrays_allocated_mem*2);
	  printf("The total ram of this system is estimated to be %zu.\n",sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGE_SIZE));
	  printf("The ram memory usable by the program according to the parameter MAX_MEMORY_ALLOC in particle_definitions.h (which sets the fraction of the available memory) is %zu\n.",max_memory_allocatable_data);
	  printf("Please, check your settings again. I quit.\n");
	  exit(4);
  } 
  avg(argv,argc,Tp,Jp,Jb,Jc,Js,Jt,Pnum,&nevents);
}
else //we already checked all main options, so we should never come here
{
  printf("Unknown option %s, strangely not detected by previous checking... Did you make, maybe, a bad code modification?\n",argv[1]);
  exit(1);
}
printf("Writing the energy momentum tensor and the four currents.\n");
write_results(outputfile,Tp,Jp,Jb,Jc,Js,Jt,Pnum,nevents);

print_dens=atoi(argv[2]);
switch(print_dens)
{
	case 0:
	printf("Not printing files with densities.\n");
	break;
	
	case 1:
	write_densities(outputfile,Tp,Jp,Jb,Jc,Js,Jt,Pnum,nevents);
	break;
	
	default:
	printf("Print option neither 0 nor 1. In the doubt, I don't compute and print the densities.\n");
}

free(time_int_array);
free(Tp);
free(Jp);
free(Jb);
free(Jc);
free(Js);
free(Jt);
free(Pnum);
#ifdef SMASH
free(plist);
#endif
big_arrays_allocated_mem=0;
return 0;
}

/** @brief help function
* It displays the program syntax.
* @param None
**/
void help()
{
printf("Syntax: ./cg <comp|avg> <0 (not print densities) | 1 (print)> <outpufile_prefix> <UrQMD_event_file_1|Tmunu_file1> [UrQMD_event_file_2|Tmunu_file2] [UrQMD_event_file_3|Tmunu_file3] ... <file_with_time_intervals|Tmunu_fileN>\n");
}





