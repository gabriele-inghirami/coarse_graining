/**
 * @file definitions.h
 *
 * @brief here we define the dimension of the grid, the data structures and we declare functions
 *
 */

#ifndef DEF_MAIN
#define DEF_MAIN

#include <inttypes.h>
#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <wchar.h>

/**
 * choice of the transport code to be used: URQMD or SMASH
 *
 */
#define URQMD
//#define SMASH

// we check that either URQMD or SMASH are defined, if not a compilation error is raised
#ifdef URQMD
#elif defined(SMASH)
#else
#error "You must define either SMASH or URQMD in definitions.h - compilation failed."
#endif

/**
 * DX_DEF the lenght (in fm) of the side x of the coarse grained cell
 *
 */
#define DX_DEF 1.0

/**
 * DY_DEF the lenght (in fm) of the side y of the coarse grained cell
 *
 */
#define DY_DEF 1.0

/**
 * DZ_DEF the lenght (in fm) of the side z of the coarse grained cell
 *
 */
#define DZ_DEF 1.0

/**
 * NX_DEF the number of coarse grained cells in the x direction
 *
 */
#define NX_DEF 7

/**
 * NY_DEF the number of coarse grained cells in the y direction
 *
 */
#define NY_DEF 7

/**
 * NZ_DEF the number of coarse grained cells in the z direction
 *
 */
#define NZ_DEF 7

/**
 * USE_CENTERED_GRID if 0 the grid is centered around 0, otherwise the grid starts at {X,Y,Z}_START
 *
 */
#define USE_CENTERED_GRID 0

/**
 * X_START the x position of the left border of first cell along the x axis
 *
 */
#define X_START -0.5

/**
 * Y_START the y position of the left border of first cell along the y axis
 *
 */
#define Y_START -0.5

/**
 * Z_START the z position of the left border of first cell along the z axis
 *
 */
#define Z_START -0.5

/**
 * B_SELECTION if not zero (e.g. 1) it selects only events with impact parameter between BMIN and BMAX, chosen below.
 * If 0 all events are accepted.
 *
 */

#define B_SELECTION 0

/**
 * BMIN the minimum impact parameter b in fm to select an event
 *
 */
#define BMIN 0

/**
 * BMAX the maximum impact parameter b in fm to select an event
 *
 */
#define BMAX 3.5

/**
 * MAX_MEMORY_ALLOC defines the maximum fraction of the total memory of the system that the program should use. If this
 * value is exceeded, the data must processed in intermediate steps.
 *
 */
#define MAX_MEMORY_ALLOC 0.98

// in a previous version we used long int for all variable to avoid type casting, we should check if there is a
// significant difference in the speed of execution
#define TLOC 10 * (p + np * (k + nz * (j + ny * (i + (long)h * nx))))
#define JPL 4 * (p + np * (k + nz * (j + ny * (i + (long)h * nx))))
#define JBL 4 * (k + nz * (j + ny * (i + (long)h * nx)))
#define PNLOC (p + np * (k + nz * (j + ny * (i + (long)h * nx))))

/**
 * NP the number of particles: from 0 to all "stable" (lifetime > 10 fm) particles (maximum 35)
 * set NP=0, leaving just the "catchall" entry, if not interested in computing the individual hadron currents
 *
 */
#define NP 35
#if NP > 35
#error NP cannot be larger than 35!
#endif

#ifdef URQMD
/**
 * \struct pdata (UrQMD)
 * It contains most of the information provided by UrQMD about a particle. However, please, note that the first member
 * is the time index instead of the actual time. The pdata components are:
 *  - t_intex : time index
 *  - x       : the x coordinate
 *  - y       : the y coordinate
 *  - z       : the z coordinate
 *  - en      : the energy
 *  - px      : the x component of the four momentum
 *  - py      : the y component of the four momentum
 *  - pz      : the z component of the four momentum
 *  - m       : the mass
 *  - itype   : the UrQMD's itype
 *  - iso3    : isospin
 *  - charge  : the electric charge
 *  - *next   : a pointer to the memory location of the next pdata structure
 */

typedef struct pdata
{
  int t_index;
  double x;
  double y;
  double z;
  double en;
  double px;
  double py;
  double pz;
  double m;
  int itype;
  int iso3;
  int charge;
  struct pdata *next;
} pdata;

#elif defined(SMASH)
/**
 * \struct pdata (SMASH)
 * It contains most of the information provided by SMASH about a particle. However, please, note that the first member
 * is the time index instead of the actual time. The pdata components are:
 *  - t_intex : time index
 *  - x       : the x coordinate
 *  - y       : the y coordinate
 *  - z       : the z coordinate
 *  - en      : the energy
 *  - px      : the x component of the four momentum
 *  - py      : the y component of the four momentum
 *  - pz      : the z component of the four momentum
 *  - m       : the mass
 *  - pdg_id  : the PDG id
 *  - charge  : the electric charge
 *  - *next   : a pointer to the memory location of the next pdata structure
 */
typedef struct pdata
{
  int t_index;
  double x;
  double y;
  double z;
  double en;
  double px;
  double py;
  double pz;
  double m;
  int pdg_id;
  int charge;
  struct pdata *next;
} pdata;
#else
#error "You must define either SMASH or URQMD in definitions.h - compilation failed."
#endif

/**
 * \struct pinfo (SMASH)
 * It contains the data of certain particle species.
 * The pinfo components are:
 *  - pdg_id     : the PDG id
 *  - name       : the name of the particle
 *  - mass       : the mass
 *  - spin       : the spin multiplicity 2J+1, J being the total spin
 *  - strangeness: the strangeness
 *  - baryon_num : the baryon number
 */
typedef struct pinfo
{
  int pdg_id;
  wchar_t name[16];
  double mass;
  int spin;
  int strangeness;
  int baryon_num;
} pinfo;

/**
 * MAX_SMASH_HADRON_SPECIES
 *
 * the maximum number of hadrons species that is expected to find in SMASH's particles.txt
 *
 */
#define MAX_SMASH_HADRON_SPECIES 500

/**
 * prepare_smash_hadron_list (SMASH)
 * It reads the file particles.txt, that comes with SMASH, and it prepares a sorted array (ordered by growind pdg_id)
 * with the properties of the hadrons managed by SMASH
 */
void prepare_smash_hadron_array ();

/** @brief (SMASH) it returns the strangeness and the baryon number of a hadron, given the pdg_id
 *
 * @param[in] the pdg_id
 *
 * @param[out] the strangeness of the hadron
 *
 * @param[out] the baryon number of the hadron
 *
 * \callgraph
 */
void get_hadron_info (const int pdg_id, int *s, int *B);

enum
{
  T00 = 0,
  T01,
  T02,
  T03,
  T11,
  T12,
  T13,
  T22,
  T23,
  T33
}; /**< \enum The flattened indexes of the energy momentum tensor components */
enum
{
  J0 = 0,
  J1,
  J2,
  J3
}; /**< \enum The flattened indexes of four momentum components */

/** @brief It computes the energy momentum tensor and the four currents
 *
 *   @param[in] an array with the names of the files to be processed
 *
 *   @param[in] an integer with the number of the files to be processed
 *
 *   @param[in] *Tp a pointer to the energy momentum tensor of the particle species
 *
 *   @param[in] *Jp a pointer to the partticle species four current
 *
 *   @param[in] *Jb a pointer to the net baryon four current
 *
 *   @param[in] *Jc a pointer to the net electric four current
 *
 *   @param[in] *Js a pointer to the net strangeness four current
 *
 *   @param[in] *Jt a pointer to the total baryon four current
 *
 *   @param[in] *Pnum a pointer to the number of hadrons for each species
 *
 *   @param[out] nevents the number of events (long int)
 */
int compute (char **, int, double_t *, double_t *, double_t *, double_t *, double_t *, double_t *, long int *,
             long int *);

/** @brief It averages the results of previously computed energy momentum tensors and four currents
 *
 *   @param[in] an array with the names of the files to be processed
 *
 *   @param[in] an integer with the number of the files to be processed
 *
 *   @param[in] *Tp a pointer to the energy momentum tensor of the particle species
 *
 *   @param[in] *Jp a pointer to the partticle species four current
 *
 *   @param[in] *Jb a pointer to the net baryon four current
 *
 *   @param[in] *Jc a pointer to the net electric charge four current
 *
 *   @param[in] *Js a pointer to the net strangeness four current
 *
 *   @param[in] *Jt a pointer to the total baryon four current
 *
 *   @param[in] *Pnum a pointer to the number of hadrons for each species
 *
 *   @param[out] nevents the number of events (long int)
 *
 */
int avg (char **, int, double_t *, double_t *, double_t *, double_t *, double_t *, double_t *, long int *, long int *);

/** @brief It processes data
 */

int process_data (double_t *, double_t *, double_t *, double_t *, double_t *, double_t *, long int *, long int,
                  pdata *);

/** @brief It returns the index of the particle in Tmunu or the number of Tmunu entries if 0,0 (UrQMD) or 0 (SMASH) are
 * given as arguments The version for UrQMD takes itype and iso3, the version for SMASH takes pdg_id Unkown particles
 * go in the last Tmunu index.
 */
#ifdef URQMD
int get_particle_index (int, int); // The third argument tells if it is a particle (+1) or an antiparticle (-1)
#elif defined(SMASH)
int get_particle_index (int);
#endif

/** @brief It returns the strangeness of an hadron
 * *
 * *   @ param[in] the itype of the particle
 * *
 * *   @ param[out] the strangeness (-(number_of_strange_quarks - number_of_anti_strange_quarks)) of the hadron
 * */
int get_strangeness (int);

/** @brief it reads the datta from an UrQMD output file or a SMASH binary output file
 *
 *   @param[in] a char pointer to the name of the file
 *
 *   @param[in] *Tp a pointer to the energy momentum tensor of the particle species
 *
 *   @param[in] *Jp a pointer to the partticle species four current
 *
 *   @param[in] *Jb a pointer to the net baryon four current
 *
 *   @param[in] *Jc a pointer to the net electric charge four current
 *
 *   @param[in] *Js a pointer to the net strangeness four current
 *
 *   @param[in] *Jt a pointer to the total baryon four current
 *
 *   @param[in] *Pnum a pointer to the number of hadrons for each species
 *
 *   @param[out] nevents the number of events (long int)
 *
 *   @param[in/out] *data_time a pointer to the system evolution time at which to perform the computation
 *
 *   \callgraph
 */

void read_data (char *, double_t *, double_t *, double_t *, double_t *, double_t *, double_t *, long int *, long int *,
                double_t *);

/** @brief it takes the timesteps to analize from a text files and it stores them into an array
 *
 *   @param[in] timefile the name of the text file containing the times to analize
 *
 *   @param[out] ntimesteps the number of timesteps
 *
 *   @param[out] timesteps_array a pointer to an array with the timesteps to analize
 *
 */
void get_timesteps (char *, int *, double **);

/** @brief it write the densities into an output file. This function is enabled/disabled when the program is invoked.
 *
 *   @param[in] outputprefix the prefix of the output files
 *
 *   @param[in] *Tp a pointer to the energy momentum tensor of the particle species
 *
 *   @param[in] *Jp a pointer to the partticle species four current
 *
 *   @param[in] *Jb a pointer to the net baryon four current
 *
 *   @param[in] *Jc a pointer to the net electric charge four current
 *
 *   @param[in] *Js a pointer to the net strangeness four current
 *
 *   @param[in] *Jt a pointer to the total baryon four current
 *
 *   @param[in] *Pnum a pointer to the number of hadrons for each species
 *
 *   @param[out] nevents the number of events (long int)
 *   \callgraph
 */
void write_densities (char *, double_t *, double_t *, double_t *, double_t *, double_t *, double_t *, long int *,
                      long int);

/** @brief it write the results into an output file. These output files contain all the necessary informations for
 * further averaging.
 *
 *   @param[in] outputprefix the prefix of the output files
 *
 *   @param[in] *Tp a pointer to the energy momentum tensor of the particle species
 *
 *   @param[in] *Jp a pointer to the particle species four current
 *
 *   @param[in] *Jb a pointer to the net baryon four current
 *
 *   @param[in] *Jc a pointer to the net electric charge four current
 *
 *   @param[in] *Js a pointer to the net strangeness four current
 *
 *   @param[in] *Jt a pointer to the total baryon four current
 *
 *   @param[in] *Pnum a pointer to the number of hadrons for each species
 *
 *   @param[in] nevents the number of events (long int)
 *   \callgraph
 */
void write_results (char *, double_t *, double_t *, double_t *, double_t *, double_t *, double_t *, long int *,
                    long int);

/** @brief it checks all input files exist. If not, it stops the program execution.
 *   \callgraph
 */
void check_input_files (char **, int);

/** @brief it returns the name of a particle from its index
 */
const char *associate_particle_array_index (int);

/** @brief it associates to a hadron it baryon number, strangeness, total spin and electric charge from its PDG id
 *
 * @param[in] the PDG id string
 *
 * @param[in] the length of the PDG id string
 *
 * @param[out] the strangeness (- number of strange quarks)
 *
 * @param[out] the baryon number
 *
 * @param[out] 2J+1, J being the total angular momentum
 * \callgraph
 */
void get_had_prop (char *, int, int *, int *, int *);

/** @brief it fills the array of hadron properties from file particles.txt (SMASH)
 *
 *  @param[in] the array of pinfo to be filled (plist)
 *
 *  @param[out] the number of hadrons registered into the array
 *
 *  \callgraph
 */
int fill (pinfo *);

/** @brief it returns the index in the array of times time_int_array corresponding to the given test_time
 */
int check_test_time (double, double *, int);

void help ();

#endif
