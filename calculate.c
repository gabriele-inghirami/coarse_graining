/* Author: Gabriele Inghirami g.inghirami@gsi.de (2019-2022) - License: GPLv3 */

#include "definitions.h"

/** @file calculate.c
 *
 *   @brief this file contains the core functions of the program: compute, avg
 * and process_data
 *
 */

extern size_t pdata_totmem, max_memory_allocatable_data, big_arrays_allocated_mem;
extern int nt, np, nr;
extern const int nx, ny, nz;
extern const double dx, dy, dz;
extern const double xmin, ymin, zmin;
extern const int start_index, index_of_output_file;
extern float_t const bmin, bmax;
extern short int b_selection;
extern double *time_int_array;
extern const char Tplabel[];
extern const char infolabel[];
extern const int shift_resonance_index;
extern const int include_total_baryon;
extern const int include_resonances;
extern const int use_urqmd;

int
compute (char **files, int ninfiles, double *Tp, double *Jp, double *Jb, double *Jc, double *Js, double *Jt,
         double *Jr, double *Tr, long int *Pnum, long int *Rnum, long int *nevents)
{
  char *fout;         // name of the output file
  int k, idx, nf, it; // number of timesteps to examine, number of input files
  FILE *infile;
  size_t buffer_size = 180;
  ssize_t numchar;
  char bitbucket[180];
  char *buffer;   // buffer for getline, i.e. to examine line by line the UrQMD
                  // .f14 files
  float b;        // impact parameter for urqmd
  double b_smash; // impact parameter for smash
  int part_ts;
  double test_time;
  int keep_searching;
  pdata *pdata_current;
  pdata *pdata_new;
  pdata *pdata_start;
  *nevents = 0;
  int take_event;
  int index_for_avg_vs_infiles = 0;
  int bb1, bb2, bb3, bb4, bb5;
  double bb0, out_time;
  int time_index;
  char bb_char[4];           // bitbucket array for the magic number in SMASH header
  char empty_char;           // character at the end of an event block in SMASH output
                             // binary
  uint16_t format_version;   // format version in SMASH header, we will check that
                             // it is 4
  uint16_t format_variant;   // format variant in SMASH header, we will check that
                             // it is 0
  uint32_t smash_vprint_len; // number of characters of the SMASH version
  char smash_ver[100];       // array for the SMASH version
  int continue_flag = 0;
  // typedef struct smash_record{double t; double  x; double  y; double z;
  // double mass; double p0; double px; double py; double pz; int32_t pdg;
  // int32_t ID; int32_t charge;} smash_record;
  const size_t smash_record_size = 9 * sizeof (double) + 3 * sizeof (int32_t);
  void *buffer_smash = NULL;             // buffer to read an output block of SMASH binary output
  void *buffer_smash_tmp = NULL;         // temporary buffer for realloc
  uint32_t elements_of_buffer_smash = 0; // dimension of the buffer used to read a block of data
  pdata *pdata_at_event_start = NULL;    // pointer to the first particle at the beginning of an event
  char p_char_at_smash_block_start;      // p character before a time block in SMASH
                                         // binary output
  uint32_t n_part_lines;                 // lines with particle information in a time block in
                                         // SMASH binary output
  fpos_t begin_of_event_file_position;   // the position of the beginning of the
                                         // current event in the opened SMASH
                                         // binary output file
  uint32_t event_number, event_number_check;
  int look_at_b;
  double bmin_sm = (float_t)bmin;
  double bmax_sm = (float_t)bmax;
  int b_selection_sm = (int)b_selection;
  size_t ret_it;
  ssize_t gl_bb; // right now, just to avoid that the compiler complains while
                 // not checking the return value of getline
  int fsf_bb;    // right now, just to avoid that the compiler complains while not
                 // checking the return value of fscanf

  const int verbose_level = 0; // set to > 0 to print advancement messages

  fout = files[index_of_output_file];
  nf = ninfiles - start_index - 1; // the number of input files with UrQMD data, excluding the last one
                                   // which contains the informations about timesteps

  if (use_urqmd)
    {
      pdata_current = (pdata *)malloc (sizeof (pdata));
      pdata_start = pdata_current;
      buffer = malloc (buffer_size * sizeof (char));
      if (buffer == NULL)
        {
          printf ("buffer allocation to examine UrQMD .f14 files failed. I quit.\n");
          exit (2);
        }

      // main loop on input files
      for (idx = 0; idx < nf; idx++)
        {
          infile = fopen (files[idx + start_index],
                          "r"); // we open the input files with data
          if (verbose_level > 0)
            {
              printf ("Working on file %s\n", files[idx + start_index]);
            }
          numchar = getline (&buffer, &buffer_size, infile);
          while (numchar > 0)
            {
              if (strncmp (buffer, "UQMD", 4) == 0)
                {
                  if (verbose_level > 0)
                    {
                      printf ("New event found\n");
                    }
                  for (k = 0; k < 3; k++)
                    gl_bb = getline (&buffer, &buffer_size, infile);
                  sscanf (buffer, "%s %f", bitbucket, &b);
                  if (((b >= bmin) && (b <= bmax)) || (b_selection == 0))
                    {
                      *nevents += 1;
                      take_event = 0;
                    }
                  else
                    take_event = -1;
                  for (k = 0; k < 14; k++)
                    gl_bb = getline (&buffer, &buffer_size,
                                     infile); // we move 14 lines forward
                }
              sscanf (buffer, "%d", &part_ts);
              gl_bb = getline (&buffer, &buffer_size, infile); // we skip another line
              // we read the line and check the timestep
              gl_bb = getline (&buffer, &buffer_size, infile);
              sscanf (buffer, "%lf", &test_time);
              time_index = check_test_time (test_time, time_int_array, nt);
              if ((time_index > -1) && (take_event == 0))
                {
                  if (verbose_level > 0)
                    {
                      printf ("time is: %lf\n", test_time);
                    }
                  pdata_new = (pdata *)malloc (sizeof (pdata));
                  pdata_totmem += (size_t)sizeof (pdata);
                  pdata_current->next = pdata_new;
                  pdata_current = pdata_new;
                  pdata_current->t_index = time_index;
                  sscanf (buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d\n", &bb0, &(pdata_current->x),
                          &(pdata_current->y), &(pdata_current->z), &(pdata_current->en), &(pdata_current->px),
                          &(pdata_current->py), &(pdata_current->pz), &(pdata_current->m),
                          &(pdata_current->itype_or_pdg_id), &(pdata_current->iso3), &(pdata_current->charge), &bb2,
                          &bb3, &bb4);
                  pdata_current->next = NULL;
                  for (k = 0; k < part_ts - 1; k++)
                    {
                      pdata_new = (pdata *)malloc (sizeof (pdata));
                      pdata_totmem += (size_t)sizeof (pdata);
                      pdata_current->next = pdata_new;
                      pdata_current = pdata_new;
                      pdata_current->t_index = time_index;
                      fsf_bb = fscanf (infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d\n", &bb0,
                                       &(pdata_current->x), &(pdata_current->y), &(pdata_current->z),
                                       &(pdata_current->en), &(pdata_current->px), &(pdata_current->py),
                                       &(pdata_current->pz), &(pdata_current->m), &(pdata_current->itype_or_pdg_id),
                                       &(pdata_current->iso3), &(pdata_current->charge), &bb2, &bb3, &bb4);
                      pdata_current->next = NULL;
                    }
                  if (verbose_level > 0)
                    {
                      printf ("Data of timestep read (%d items)\n", part_ts);
                    }
                }
              else
                {
                  for (k = 0; k < part_ts - 1; k++)
                    gl_bb = getline (&buffer, &buffer_size, infile);
                }
              numchar = getline (&buffer, &buffer_size, infile);
            }
          fclose (infile);
          printf ("File %s read.\n", files[idx + start_index]);
          if (((pdata_totmem + big_arrays_allocated_mem) > max_memory_allocatable_data)
              && (idx < nf - 1)) // if we are working on the last file, we will process the data anyway
            {
              process_data (Tp, Jp, Jb, Jc, Js, Jt, Jt, Tr, Pnum, Rnum, *nevents, pdata_start);
              // process_data, hopefully, should deallocate all the memory for pdata
              // structures in the linked list, so we reset the memory counter and we
              // repeat the initial allocation step
              pdata_totmem = 0;
              pdata_current = (pdata *)malloc (sizeof (pdata));
              pdata_start = pdata_current;
            }
        }

      // we free the allocated arrays
      free (buffer);
    }
  else // SMASH case
    {
      fout = files[index_of_output_file];
      nf = ninfiles - start_index - 1; // the number of input files with UrQMD data, excluding the last one
                                       // which contains the informations about timesteps

      pdata_current = (pdata *)malloc (sizeof (pdata));
      pdata_start = pdata_current;

      // main loop on input files
      for (idx = 0; idx < nf; idx++)
        {
          infile = fopen (files[idx + start_index],
                          "rb"); // we open the input files with data
          printf ("Working on file %s\n", files[idx + start_index]);
          ret_it = ret_it = fread (bb_char, sizeof (char), 4, infile);
          if (ret_it == 0)
            {
              printf ("Data reading failure. Exiting.\n");
              exit (4);
            }
          if (ret_it == 0)
            {
              printf ("Data reading failure. Exiting.\n");
              exit (4);
            }
          ret_it = fread (&format_version, sizeof (uint16_t), 1, infile);
          if (ret_it == 0)
            {
              printf ("Data reading failure. Exiting.\n");
              exit (4);
            }
          /*
          if((int)format_version!=4)
          {
                  printf("Error! This program is designed to read SHASH binary output
          version 4, but you provided a file written according to version
          %u\n.",format_version); exit(7);
          }*/
          ret_it = fread (&format_variant, sizeof (uint16_t), 1, infile);
          if (ret_it == 0)
            {
              printf ("Data reading failure. Exiting.\n");
              exit (4);
            }
          /*if(format_variant!=0)
          {
                  printf("Error! This program is designed to read SHASH format variant
          0 (normal), but you provided a file written according to variant
          %u\n.",format_variant); exit(7);
          }*/
          ret_it = fread (&smash_vprint_len, sizeof (uint32_t), 1, infile);
          if (ret_it == 0)
            {
              printf ("Data reading failure. Exiting.\n");
              exit (4);
            }
          ret_it = fread (smash_ver, sizeof (char), smash_vprint_len, infile);
          if (ret_it == 0)
            {
              printf ("Data reading failure. Exiting.\n");
              exit (4);
            }
          smash_ver[smash_vprint_len] = '\0';
          printf ("Reading %s binary output, OSCAR format version %u, variant %u\n", smash_ver, format_version,
                  format_variant);
          fgetpos (infile, &begin_of_event_file_position);
          if (b_selection_sm == 0)
            look_at_b = 0;
          else
            look_at_b = 1;
          event_number_check = 0;
          while (1)
            {
              if (ret_it = fread (&p_char_at_smash_block_start, sizeof (char), 1, infile) < 1)
                break;
              if (ret_it == 0)
                {
                  printf ("Data reading failure. Exiting.\n");
                  exit (4);
                }
              if (!((p_char_at_smash_block_start == 'p') || (p_char_at_smash_block_start == 'f')))
                {
                  printf ("Error in reading file %s. I did not find the p character at "
                          "the beginning of a block header, "
                          "but %c. I quit.\n",
                          files[idx + start_index], p_char_at_smash_block_start);
                  exit (8);
                }
              if (p_char_at_smash_block_start == 'p')
                {
                  ret_it = fread (&n_part_lines, sizeof (uint32_t), 1, infile);
                  if (ret_it == 0)
                    {
                      printf ("Data reading failure. Exiting.\n");
                      exit (4);
                    }
                  if ((b_selection_sm != 0) && (look_at_b == 1))
                    {
                      // we skip the lines, we want to go a the end of the file and read the
                      // value of the impact parameter
                      fseek (infile, n_part_lines * smash_record_size, SEEK_CUR);
                    }
                  else
                    {
                      // we read the first double, to check the time
                      ret_it = fread (&out_time, sizeof (double), 1, infile);
                      if (ret_it == 0)
                        {
                          printf ("Data reading failure. Exiting.\n");
                          exit (4);
                        }
                      fseek (infile, -sizeof (double), SEEK_CUR); // we move back
                      time_index = check_test_time (out_time, time_int_array, nt);
                      if (time_index < 0) // the time in this block is not among those to be analized
                        {
                          fseek (infile, n_part_lines * smash_record_size, SEEK_CUR);
                        }
                      else
                        {
                          // we read a block of particles at time, we check if the impact
                          // parameter is within the desired range and, if so, we analize the
                          // data but first, we check if we have allocated enough room
                          if (buffer_smash == NULL)
                            {
                              buffer_smash = malloc (smash_record_size * n_part_lines);
                              elements_of_buffer_smash = n_part_lines;
                              if (buffer_smash == NULL)
                                {
                                  printf ("Unable to allocate buffer_smash for output block of "
                                          "file %s. I quit.\n",
                                          files[idx + start_index]);
                                  exit (8);
                                }
                            }
                          else if (n_part_lines > elements_of_buffer_smash)
                            {
                              buffer_smash_tmp = realloc ((void *)buffer_smash, smash_record_size * n_part_lines);
                              if (buffer_smash == NULL)
                                {
                                  printf ("Unable to reallocate buffer_smash for output block of "
                                          "file %s. I quit.\n",
                                          files[idx + start_index]);
                                  exit (8);
                                }
                              buffer_smash = buffer_smash_tmp;
                              buffer_smash_tmp = NULL;
                              elements_of_buffer_smash = n_part_lines;
                            }
                          ret_it = fread (buffer_smash, smash_record_size, n_part_lines, infile);
                          if (ret_it == 0)
                            {
                              printf ("Data reading failure. Exiting.\n");
                              exit (4);
                            }
                          // we have read the data, now we store them into the pdata
                          // linked list
                          for (k = 0; k < n_part_lines; k++)
                            {
                              pdata_new = (pdata *)malloc (sizeof (pdata));
                              pdata_totmem += (size_t)sizeof (pdata);
                              pdata_current->next = pdata_new;
                              pdata_current = pdata_new;
                              pdata_current->t_index = time_index;
                              // we read the data from the buffer
                              buffer_smash += sizeof (double); // we skip the first information about time
                              memcpy (&(pdata_current->x), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->y), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->z), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->m), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->en), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->px), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->py), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy (&(pdata_current->pz), buffer_smash, sizeof (double));
                              buffer_smash += sizeof (double); // we move forward 1 double
                              memcpy ((int *)&(pdata_current->itype_or_pdg_id), buffer_smash, sizeof (uint32_t));
                              buffer_smash += 2 * sizeof (uint32_t); // we move fwd 2 unsigned int, because
                                                                     // we skip the particle internal ID info
                              memcpy ((int *)&(pdata_current->charge), buffer_smash, sizeof (uint32_t));
                              buffer_smash += sizeof (uint32_t); // we move fwd 1 unsigned int

                              pdata_current->next = NULL;
                            }
                          // printf("Read %d lines \n",n_part_lines);
                          // now we restore buffer_smash to its original position
                          buffer_smash = buffer_smash_tmp;
                          buffer_smash_tmp = NULL;
                        }
                    }
                }
              else // clearly, p_char_at_smash_block_start is now 'f'
                {
                  ret_it = fread (&event_number, sizeof (uint32_t), 1, infile);
                  if (ret_it == 0)
                    {
                      printf ("Data reading failure. Exiting.\n");
                      exit (4);
                    }
                  if (event_number != event_number_check)
                    {
                      printf ("There is a discrepancy between the number of events "
                              "according to the file %s (%u) and what "
                              "we counted so far (%u)... Please, check!\n",
                              files[idx + start_index], event_number, event_number_check);
                      exit (8);
                    }
                  else
                    {
                      if (look_at_b == 0)
                        {
                          event_number_check = event_number_check + 1;
                        }
                    }
                  ret_it = fread (&b_smash, sizeof (double), 1, infile);
                  if (ret_it == 0)
                    {
                      printf ("Data reading failure. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&empty_char, sizeof (char), 1, infile);
                  if (ret_it == 0)
                    {
                      printf ("Data reading failure. Exiting.\n");
                      exit (4);
                    }
                  if (((b_smash >= bmin_sm) && (b_smash <= bmax_sm)) || (b_selection_sm == 0))
                    {
                      if (look_at_b == 0)
                        {
                          *nevents = *nevents + 1;
                          if (b_selection_sm != 0)
                            {
                              look_at_b = 1;
                              fgetpos (infile, &begin_of_event_file_position);
                            }
                        }
                      else
                        {
                          if (b_selection_sm != 0)
                            {
                              look_at_b = 0;
                              fsetpos (infile, &begin_of_event_file_position);
                            }
                        }
                    }
                  else
                    {
                      event_number_check = event_number_check + 1; // we discard the event, but we update the
                                                                   // event_number_check variable
                      fgetpos (infile, &begin_of_event_file_position);
                    }
                }
            }
          fclose (infile);
          printf ("File %s read, with events %u.\n", files[idx + start_index], event_number + 1);
          if (((pdata_totmem + big_arrays_allocated_mem) > max_memory_allocatable_data)
              && (idx < nf - 1)) // if we are working on the last file, we will process the data anyway
            {
              process_data (Tp, Jp, Jb, Jc, Js, Jt, Jr, Tr, Pnum, Rnum, *nevents, pdata_start);
              // process_data, hopefully, should deallocate all the memory for pdata
              // structures in the linked list, so we reset the memory counter and we
              // repeat the initial allocation step
              pdata_totmem = 0;
              pdata_current = (pdata *)malloc (sizeof (pdata));
              pdata_start = pdata_current;
            }
        }

      free (buffer_smash);
      buffer_smash = NULL;
    }

  process_data (Tp, Jp, Jb, Jc, Js, Jt, Jr, Tr, Pnum, Rnum, *nevents, pdata_start);

  printf ("Number of events: %ld\n", *nevents);

  return 0;
}

int
avg (char **files, int ninfiles, double *Tp, double *Jp, double *Jb, double *Jc, double *Js, double *Jt, double *Jr,
     double *Tr, long int *Pnum, long int *Rnum, long int *nevents)
{
  double *Tp2, *Jp2, *Jb2, *Jc2, *Js2, *Jt2, *Tr2, *Jr2;
  long int *Pnum2, *Rnum2;
  double time_check;
  int nf, h, f, i, j, k, l, p, r;
  long int nevents2;
  char *input_filename_Tp;
  char time_suffix[8];
  int prefixlen;
  char *inputfile;
  double timeref;

  timeref = -1;
  h = 0; // h is used in the TLOC expression, defined in definitions.h

  Tp2 = (double *)malloc (sizeof (double) * nx * ny * nz * np * 10);
  if (Tp2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Tp2 array inside "
              "avg. I am forced to quit.\n");
      exit (4);
    }
  Jp2 = (double *)malloc (sizeof (double) * nx * ny * nz * np * 4);
  if (Jp2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jp2 array inside "
              "avg. I am forced to quit.\n");
      exit (4);
    }
  Jb2 = (double *)malloc (sizeof (double) * nx * ny * nz * 4);
  if (Jb2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jb2 array inside "
              "main. I am forced to quit.\n");
      exit (4);
    }
  Jc2 = (double *)malloc (sizeof (double) * nx * ny * nz * 4);
  if (Jc2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jc2 array inside "
              "main. I am forced to quit.\n");
      exit (4);
    }
  Js2 = (double *)malloc (sizeof (double) * nx * ny * nz * 4);
  if (Js2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Js2 array inside "
              "main. I am forced to quit.\n");
      exit (4);
    }
  if (include_total_baryon)
    {
      Jt2 = (double *)malloc (sizeof (double) * nx * ny * nz * 4);
      if (Jt2 == NULL)
        {
          printf ("Sorry, but it is not possible to allocate the Jt2 array inside "
                  "main. I am forced to quit.\n");
          exit (4);
        }
    }
  else
    {
      Jt2 = (double *)malloc (sizeof (double));
    }
  if (include_resonances)
    {
      Jr2 = (double *)malloc (sizeof (double) * nx * ny * nz * nr * 4);
      if (Jr2 == NULL)
        {
          printf ("Sorry, but it is not possible to allocate the Jr2 array inside "
                  "avg. I am forced to quit.\n");
          exit (4);
        }
      Tr2 = (double *)malloc (sizeof (double) * nx * ny * nz * nr * 10);
      if (Tr2 == NULL)
        {
          printf ("Sorry, but it is not possible to allocate the Tr2 array inside "
                  "avg. I am forced to quit.\n");
          exit (4);
        }
      Rnum2 = (long int *)malloc (sizeof (long int) * nx * ny * nz * nr);
      if (Rnum2 == NULL)
        {
          printf ("Sorry, but it is not possible to allocate the Rnum2 array inside "
                  "main. I am forced to quit.\n");
          exit (4);
        }
    }
  else
    {
      Jr2 = (double *)malloc (sizeof (double));
      Tr2 = (double *)malloc (sizeof (double));
      Rnum2 = (long int *)malloc (sizeof (long int));
    }

  Pnum2 = (long int *)malloc (sizeof (long int) * nx * ny * nz * np);
  if (Pnum2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Pnum2 array inside "
              "main. I am forced to quit.\n");
      exit (4);
    }
  // we already checked that the allocated memory was less than the maximum
  // allocatable in the first part of the main function, so we do not do it
  // again here
  nf = ninfiles - start_index; // the number of input files with Tmunu data

  for (f = 0; f < nf; f++)
    {
      inputfile = files[start_index + f];

      read_data (inputfile, Tp2, Jp2, Jb2, Jc2, Js2, Jt2, Jr2, Tr2, Pnum2, Rnum2, &nevents2, &time_check);
      if (timeref < 0)
        {
          timeref = time_check;
          time_int_array[0] = timeref;
        }
      else
        {
          if (time_check != timeref)
            {
              printf ("Mismatching between the time in file %s (%lf) and "
                      "corresponding entry %d in time array "
                      "obtained from file %s.\n",
                      input_filename_Tp, time_check, h, files[ninfiles - 1]);
              printf ("I quit.\n");
              exit (5);
            }
        }
      *nevents += nevents2;
      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {
                  for (p = 0; p < np; p++)
                    {
                      for (r = 0; r < 10; r++)
                        Tp[r + TLOC] += Tp2[r + TLOC];
                    }
                }
            }
        }
      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {
                  for (p = 0; p < np; p++)
                    {
                      for (r = 0; r < 4; r++)
                        Jp[r + JPL] += Jp2[r + JPL];
                    }
                }
            }
        }
      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {
                  for (r = 0; r < 4; r++)
                    Jb[r + JBL] += Jb2[r + JBL];
                }
            }
        }
      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {
                  for (r = 0; r < 4; r++)
                    Jc[r + JBL] += Jc2[r + JBL];
                }
            }
        }
      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {
                  for (r = 0; r < 4; r++)
                    Js[r + JBL] += Js2[r + JBL];
                }
            }
        }
      if (include_total_baryon)
        {
          for (i = 0; i < nx; i++)
            {
              for (j = 0; j < ny; j++)
                {
                  for (k = 0; k < nz; k++)
                    {
                      for (r = 0; r < 4; r++)
                        Jt[r + JBL] += Jt2[r + JBL];
                    }
                }
            }
        }
      if (include_resonances)
        {
          for (i = 0; i < nx; i++)
            {
              for (j = 0; j < ny; j++)
                {
                  for (k = 0; k < nz; k++)
                    {
                      for (r = 0; r < nr; r++)
                        {
                          for (l = 0; l < 4; l++)
                            Jr[l + JRL] += Jr2[l + JRL];
                        }
                    }
                }
            }
          for (i = 0; i < nx; i++)
            {
              for (j = 0; j < ny; j++)
                {
                  for (k = 0; k < nz; k++)
                    {
                      for (r = 0; r < nr; r++)
                        {
                          for (l = 0; l < 10; l++)
                            Tr[l + RTLOC] += Tr2[l + RTLOC];
                        }
                    }
                }
            }
          for (i = 0; i < nx; i++)
            {
              for (j = 0; j < ny; j++)
                {
                  for (k = 0; k < nz; k++)
                    {
                      for (r = 0; r < nr; r++)
                        {
                          Rnum[RNLOC] += Rnum2[RNLOC];
                        }
                    }
                }
            }
        }
      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {
                  for (p = 0; p < np; p++)
                    {
                      Pnum[PNLOC] += Pnum2[PNLOC];
                    }
                }
            }
        }
    }
  free (Tp2);
  free (Jp2);
  free (Jb2);
  free (Jc2);
  free (Js2);
  free (Jt2);
  free (Jr2);
  free (Tr2);
  free (Pnum2);
  free (Rnum2);
  return 0;
}

int
process_data (double *Tp, double *Jp, double *Jb, double *Jc, double *Js, double *Jt, double *Jr, double *Tr,
              long int *Pnum, long int *Rnum, long int nevents, pdata *pdata_input)
{
  int h, i, j, k, p, r;
  int jsign;
  int continue_flag;
  int strangeness, baryon_number;

  /* for debugging purposes, tot1 total number of hadrons, tot2 hadrons within
  the coarse grained grid int tot1, tot2; tot1=0; tot2=0;
  */

  const int verbose_level = 0;

  pdata *pdata_entry, *pdata_next;
  pdata_entry = pdata_input->next;
  free (pdata_input);
  if (pdata_entry != NULL)
    {
      continue_flag = 0;
    }
  else
    {
      printf ("I just start computing the energy-momentum tensor, but I have no "
              "particle data. Something went "
              "wrong... I quit.\n");
      exit (2);
    }
  while (continue_flag == 0)
    {
      // tot1=tot1+1;
      i = (int)floor ((pdata_entry->x - xmin) / dx);
      if (!((i >= 0) && (i < nx)))
        {

          if (pdata_entry->next != NULL)
            {
              pdata_next = pdata_entry->next;
              free (pdata_entry);
              pdata_entry = pdata_next;
            }
          else
            {
              free (pdata_entry);
              continue_flag = 1;
            }
          continue;
        }
      j = (int)floor ((pdata_entry->y - ymin) / dy);
      if (!((j >= 0) && (j < ny)))
        {
          if (pdata_entry->next != NULL)
            {
              pdata_next = pdata_entry->next;
              free (pdata_entry);
              pdata_entry = pdata_next;
            }
          else
            {
              free (pdata_entry);
              continue_flag = 1;
            }
          continue;
        }
      k = (int)floor ((pdata_entry->z - zmin) / dz);
      if (!((k >= 0) && (k < nz)))
        {
          if (pdata_entry->next != NULL)
            {
              pdata_next = pdata_entry->next;
              free (pdata_entry);
              pdata_entry = pdata_next;
            }
          else
            {
              free (pdata_entry);
              continue_flag = 1;
            }
          continue;
        }
      h = pdata_entry->t_index;
      if (np == 0)
        {
          p = 0;
        }
      else
        {
          p = get_particle_index (pdata_entry->itype_or_pdg_id, pdata_entry->iso3);
        }
      if (include_resonances)
        {
          if (p >= shift_resonance_index)
            {
              r = p - shift_resonance_index; // we remove the offset
              p = np - 1;                    // we set the resonance index to the catch-all entry
              Rnum[RNLOC] += 1;
              Tr[T00 + RTLOC] += pdata_entry->en;
              Tr[T01 + RTLOC] += pdata_entry->px;
              Tr[T02 + RTLOC] += pdata_entry->py;
              Tr[T03 + RTLOC] += pdata_entry->pz;
              Tr[T11 + RTLOC] += (pdata_entry->px * pdata_entry->px) / (pdata_entry->en);
              Tr[T12 + RTLOC] += (pdata_entry->px * pdata_entry->py) / (pdata_entry->en);
              Tr[T13 + RTLOC] += (pdata_entry->px * pdata_entry->pz) / (pdata_entry->en);
              Tr[T22 + RTLOC] += (pdata_entry->py * pdata_entry->py) / (pdata_entry->en);
              Tr[T23 + RTLOC] += (pdata_entry->py * pdata_entry->pz) / (pdata_entry->en);
              Tr[T33 + RTLOC] += (pdata_entry->pz * pdata_entry->pz) / (pdata_entry->en);
              Jr[J0 + JRL] += 1;
              Jr[J1 + JRL] += (pdata_entry->px) / pdata_entry->en;
              Jr[J2 + JRL] += (pdata_entry->py) / pdata_entry->en;
              Jr[J3 + JRL] += (pdata_entry->pz) / pdata_entry->en;
            }
        }
      // tot2=tot2+1;
      Pnum[PNLOC] += 1;
      Tp[T00 + TLOC] += pdata_entry->en;
      Tp[T01 + TLOC] += pdata_entry->px;
      Tp[T02 + TLOC] += pdata_entry->py;
      Tp[T03 + TLOC] += pdata_entry->pz;
      Tp[T11 + TLOC] += (pdata_entry->px * pdata_entry->px) / (pdata_entry->en);
      Tp[T12 + TLOC] += (pdata_entry->px * pdata_entry->py) / (pdata_entry->en);
      Tp[T13 + TLOC] += (pdata_entry->px * pdata_entry->pz) / (pdata_entry->en);
      Tp[T22 + TLOC] += (pdata_entry->py * pdata_entry->py) / (pdata_entry->en);
      Tp[T23 + TLOC] += (pdata_entry->py * pdata_entry->pz) / (pdata_entry->en);
      Tp[T33 + TLOC] += (pdata_entry->pz * pdata_entry->pz) / (pdata_entry->en);
      Jp[J0 + JPL] += 1;
      Jp[J1 + JPL] += (pdata_entry->px) / pdata_entry->en;
      Jp[J2 + JPL] += (pdata_entry->py) / pdata_entry->en;
      Jp[J3 + JPL] += (pdata_entry->pz) / pdata_entry->en;
      if (use_urqmd)
        {
          // printf("Here\n");
          strangeness = get_strangeness (pdata_entry->itype_or_pdg_id);
          if (abs (pdata_entry->itype_or_pdg_id) < 100)
            {
              if (pdata_entry->itype_or_pdg_id > 0)
                {
                  jsign = 1;
                }
              else
                {
                  jsign = -1;
                }
            }
          else
            {
              jsign = 0;
            }
        }
      else
        {
          get_hadron_info (pdata_entry->itype_or_pdg_id, &strangeness, &baryon_number);
          // printf("hadron: %d p: %d strangeness: %d B:
          // %d\n",pdata_entry->itype_or_pdg_id,p,strangeness, baryon_number);
          if (baryon_number != 0)
            {
              jsign = baryon_number;
            }
          else
            {
              jsign = 0;
            }
        }
      if (jsign != 0)
        {
          // printf("**B %lf  %lf  %lf  %lf\n",Jb[J0 + JBL], Jb[J1 + JBL], Jb[J2 + JBL], Jb[J3 + JBL]);
          // printf("jsign: %d\n",jsign);
          Jb[J0 + JBL] += jsign;
          Jb[J1 + JBL] += jsign * (pdata_entry->px) / pdata_entry->en;
          Jb[J2 + JBL] += jsign * (pdata_entry->py) / pdata_entry->en;
          Jb[J3 + JBL] += jsign * (pdata_entry->pz) / pdata_entry->en;
          // printf("**A %lf  %lf  %lf  %lf\n",Jb[J0 + JBL], Jb[J1 + JBL], Jb[J2 + JBL], Jb[J3 + JBL]);
          if (include_total_baryon)
            {
              Jt[J0 + JBL] += 1;
              Jt[J1 + JBL] += (pdata_entry->px) / pdata_entry->en;
              Jt[J2 + JBL] += (pdata_entry->py) / pdata_entry->en;
              Jt[J3 + JBL] += (pdata_entry->pz) / pdata_entry->en;
            }
        }
      Jc[J0 + JBL] += (pdata_entry->charge);
      Jc[J1 + JBL] += (pdata_entry->charge) * (pdata_entry->px) / pdata_entry->en;
      Jc[J2 + JBL] += (pdata_entry->charge) * (pdata_entry->py) / pdata_entry->en;
      Jc[J3 + JBL] += (pdata_entry->charge) * (pdata_entry->pz) / pdata_entry->en;
      // printf("hadron: %d p: %d strangeness: %d B:
      // %d\n",pdata_entry->itype_or_pdg_id,p,strangeness, baryon_number);
      Js[J0 + JBL] += strangeness;
      Js[J1 + JBL] += strangeness * (pdata_entry->px) / pdata_entry->en;
      Js[J2 + JBL] += strangeness * (pdata_entry->py) / pdata_entry->en;
      Js[J3 + JBL] += strangeness * (pdata_entry->pz) / pdata_entry->en;

      if (pdata_entry->next != NULL)
        {
          pdata_next = pdata_entry->next;
          free (pdata_entry);
          pdata_entry = pdata_next;
        }
      else
        {
          free (pdata_entry);
          continue_flag = 1;
        }
      continue;
    }
  if (verbose_level > 0)
    {
      printf ("Four current and energy momentum tensor arrays computed\n");
    }
  // printf("tot1: %d, tot2: %d\n",tot1,tot2);
  return 0;
}
