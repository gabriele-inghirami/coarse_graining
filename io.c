/* Author: Gabriele Inghirami g.inghirami@gsi.de (2019-2022) - License: GPLv3 */

#include "definitions.h"

/** @file io.c
 *
 *   @brief this file contains the functions which write from/to disk:
 * write_results, write_densities and read_data
 */

extern int nt, np, nr;
extern const int nx, ny, nz;
extern const double dx, dy, dz;
extern const double xmin, ymin, zmin;
extern const char Tplabel[];
extern const char infolabel[];
extern const char dens_label[];
extern const char info_dens_label[];
extern double *time_int_array;
extern const int T10, T20, T30, T21, T31, T32;
extern const double cell_volume;
extern const int output_content_info;
extern const int include_total_baryon;
extern const int include_resonances;
extern const int use_urqmd;

size_t ret_it;

void
write_results (char *outputprefix, double *Tp, double *Jp, double *Jb, double *Jc, double *Js, double *Jt, double *Jr,
               double *Tr, long int *Pnum, long int *Rnum, long int nevents)
{
  int h, i, j, k, p, r;
  char *output_filename_Tp;
  char *info_filename;
  char time_suffix[8];
  static int prefixlen, Tmunu_label_len, info_label_len;
  FILE *fTp, *finfo;
  int static first_time = 0; // the first time we print also the legend and the grid spacing

  if (first_time == 0) // we print the grid informations
    {
      prefixlen = strlen (outputprefix);
      Tmunu_label_len = strlen (Tplabel);
      info_label_len = strlen (infolabel);
      info_filename
          = (char *)malloc (sizeof (char) * (prefixlen + info_label_len + 1)); // we need to add the \0 character
      *info_filename = '\0';
      strcat (info_filename, outputprefix);
      strcat (info_filename, infolabel);
      finfo = fopen (info_filename, "w");
      fprintf (finfo, "Grid order - row major order (C):\n x_0 y_0 z_0\n x_0 y_0 "
                      "z_1\n ...\n x_0 y_0 z_n-1\n x_0 y_1 "
                      "z_0\n ...\n\n");
      fprintf (finfo, "Index:   X coordinate:\n");
      for (i = 0; i < nx; i++)
        fprintf (finfo, " %4d       %9.4lf\n", i, xmin + (i + 0.5) * dx);
      fprintf (finfo, "\n");
      fprintf (finfo, "Index:   Y coordinate:\n");
      for (j = 0; j < ny; j++)
        fprintf (finfo, " %4d       %9.4lf\n", j, ymin + (j + 0.5) * dy);
      fprintf (finfo, "\n");
      fprintf (finfo, "Index:   Z coordinate:\n");
      for (k = 0; k < nz; k++)
        fprintf (finfo, " %4d       %9.4lf\n", k, zmin + (k + 0.5) * dz);
      fprintf (finfo, "\n");
      fprintf (finfo,
               "Volume factor (all saved quantities should be divided by this "
               "number): %lf\n\n",
               dx * dy * dz);
      fprintf (finfo, "Indexes for Tmunu of hadrons (%d entries):\n", np);
      for (i = 0; i < np - 1; i++)
        {
          fprintf (finfo, "%4d:\t%s\n", i, associate_particle_array_index (i));
        }
      fprintf (finfo, "%4d:\tAll other particles\n", i);
      if (include_resonances)
        {
          fprintf (finfo, "Indexes for Tmunu of resonances (%d entries):\n", nr);
          for (i = 0; i < nr; i++)
            {
              fprintf (finfo, "%4d:\t%s\n", i, associate_resonance_array_index (i));
            }
        }
      fprintf (finfo, "\nContents of the output binary file: \n");
      fprintf (finfo, "Events (long int):    %ld\n", nevents);
      fprintf (finfo, "Time (double)\n");
      fprintf (finfo, "Output content flag:    %d\n", output_content_info);
      fprintf (finfo, "Number of hadron species np (int):    %d\n", np);
      fprintf (finfo, "x cells nx (int):    %d\n", nx);
      fprintf (finfo, "y cells ny (int):    %d\n", ny);
      fprintf (finfo, "z cells nz (int):    %d\n", nz);
      fprintf (finfo, "dx (double):    %lf\n", dx);
      fprintf (finfo, "dy (double):    %lf\n", dy);
      fprintf (finfo, "dz (double):    %lf\n", dz);
      fprintf (finfo, "xmin, i.e. low border of the 1st x cell (double):    %lf\n", xmin);
      fprintf (finfo, "ymin, i.e. low border of the 1st y cell (double):    %lf\n", ymin);
      fprintf (finfo, "zmin, i.e. low border of the 1st z cell (double):    %lf\n", zmin);
      fprintf (finfo, "output_content_info (int), optional output flag:    %d\n", output_content_info);
      fprintf (finfo, "Tp, hadron enery momentum tensor, type double, size nx*ny*nz*np*10\n");
      fprintf (finfo, "Jp, hadron four current, type double, size nx*ny*nz*np*4\n");
      fprintf (finfo, "Jb, net baryon four current, type double, size nx*ny*nz*np*4\n");
      fprintf (finfo, "Jc, electric charge four current, type double, size nx*ny*nz*np*4\n");
      fprintf (finfo, "Js, strangeness four current, type double, size nx*ny*nz*np*4\n");
      if (include_total_baryon)
        {
          fprintf (finfo, "Jt, total baryon four current, type double, size nx*ny*nz*np*4\n");
        }
      fprintf (finfo, "Jt placeholder, type double, size 1\n");
      if (include_resonances)
        {
          fprintf (finfo, "Number of resonance species nr (int):    %d\n", nr);
          fprintf (finfo, "Jr, resonance four current, type double, size nx*ny*nz*np*4\n");
          fprintf (finfo, "Tr, resonance enery momentum tensor, type double, size "
                          "nx*ny*nz*nr*10\n");
        }
      else
        {
          fprintf (finfo, "Jr placeholder, type double, size 1\n");
          fprintf (finfo, "Tr placeholder, type double, size 1\n");
        }
      fprintf (finfo, "Pnum, hadron count, type double, size nx*ny*nz*np\n");
      if (include_resonances)
        {
          fprintf (finfo, "Rnum, resonance count, type double, size nx*ny*nz*np\n");
        }
      else
        {
          fprintf (finfo, "Rnum placeholder, type double, size 1\n");
        }

      fclose (finfo);
      printf ("Written grid and particle informations in file %s\n", info_filename);
      first_time = 1;
      free (info_filename);
    }

  for (h = 0; h < nt; h++)
    {
      snprintf (time_suffix, 8, "%07.3lf", time_int_array[h]);
      strncpy (&time_suffix[3], "_", 1); // we replace the dot with an underscore
      output_filename_Tp
          = (char *)malloc (sizeof (char) * (prefixlen + Tmunu_label_len + 8 + 1)); // we need to add the \0 character
      *output_filename_Tp = '\0';
      strncat (output_filename_Tp, outputprefix, prefixlen);
      strncat (output_filename_Tp, Tplabel, Tmunu_label_len);
      strncat (output_filename_Tp, time_suffix, 8);

      fTp = fopen (output_filename_Tp, "wb");
      if (fTp == NULL)
        {
          printf ("Sorry, I could not open the output file %s, so I am forced to "
                  "quit.\n",
                  output_filename_Tp);
          exit (3);
        }
      fwrite (&nevents, sizeof (long int), 1, fTp);
      fwrite (&time_int_array[h], sizeof (double), 1, fTp);
      fwrite (&output_content_info, sizeof (int), 1, fTp);
      fwrite (&np, sizeof (int), 1, fTp);
      fwrite (&nx, sizeof (int), 1, fTp);
      fwrite (&ny, sizeof (int), 1, fTp);
      fwrite (&nz, sizeof (int), 1, fTp);
      fwrite (&dx, sizeof (double), 1, fTp);
      fwrite (&dy, sizeof (double), 1, fTp);
      fwrite (&dz, sizeof (double), 1, fTp);
      fwrite (&xmin, sizeof (double), 1, fTp);
      fwrite (&ymin, sizeof (double), 1, fTp);
      fwrite (&zmin, sizeof (double), 1, fTp);
      fwrite (&Tp[h * nx * ny * nz * np * 10], sizeof (double), nx * ny * nz * np * 10, fTp);
      fwrite (&Jp[h * nx * ny * nz * np * 4], sizeof (double), nx * ny * nz * np * 4, fTp);
      fwrite (&Jb[h * nx * ny * nz * 4], sizeof (double), nx * ny * nz * 4, fTp);
      fwrite (&Jc[h * nx * ny * nz * 4], sizeof (double), nx * ny * nz * 4, fTp);
      fwrite (&Js[h * nx * ny * nz * 4], sizeof (double), nx * ny * nz * 4, fTp);
      if (include_total_baryon)
        {
          fwrite (&Jt[h * nx * ny * nz * 4], sizeof (double), nx * ny * nz * 4, fTp);
        }
      if (include_resonances)
        {
          fwrite (&nr, sizeof (int), 1, fTp);
          fwrite (&Jr[h * nx * ny * nz * nr * 4], sizeof (double), nx * ny * nz * nr * 4, fTp);
          fwrite (&Tr[h * nx * ny * nz * nr * 10], sizeof (double), nx * ny * nz * nr * 10, fTp);
        }
      fwrite (&Pnum[h * nx * ny * nz * np], sizeof (long int), nx * ny * nz * np, fTp);
      if (include_resonances)
        {
          fwrite (&Rnum[h * nx * ny * nz * nr], sizeof (long int), nx * ny * nz * nr, fTp);
        }
      printf ("Output data saved in file %s.\n", output_filename_Tp);
      fflush (fTp);
      fclose (fTp);
      free (output_filename_Tp);
    }
}

void
write_densities (char *outputprefix, double *Tp, double *Jp, double *Jb, double *Jc, double *Js, double *Jt,
                 double *Jr, double *Tr, long int *Pnum, long int *Rnum, long int nevents)
{
  int h, i, j, k, p, l, r;
  char *output_filename_Tp;
  char *info_filename;
  char time_suffix[8];
  static int prefixlen, dens_label_len, info_label_len;
  FILE *fTp, *finfo;
  int static first_time = 0; // the first time we print also the legend and the grid spacing
  char string_particle_description[2048];
  double tmp_value;
  double *empty_array; // array filled with zeroes
  /*
    the decomposition is: J^\mu (CF) = rho u^\mu (CF) + I_diff^\mu,
    where I_diff^\mu = Delta^\mu_nu J^\nu = ( Kron_delta^\mu_\nu - u^\mu u_\nu )
    J^\nu
  */
  double rho_c, rho_s, rho_t;
  double vel_B[3], vel_c[3], vel_s[3], vel_h[3];
  const double vel_null[3] = { 0, 0, 0 };
  const double null_2[2] = { 0, 0 };
  double u4[4];
  double Ic_diffusion[4], Is_diffusion[4];
  double Ic_diffcheck[4], Is_diffcheck[4];
  double jf, jf_arg, jh_arg, rho, eps;
  int offset = 14;

// uncomment the following line to print some debugging messages when printing densities
#define DBG_DENS

  if (first_time == 0) // we print the grid informations
    {
      prefixlen = strlen (outputprefix);
      dens_label_len = strlen (dens_label);
      info_label_len = strlen (info_dens_label);
      info_filename = (char *)malloc (sizeof (char) * (prefixlen + info_label_len + 1));
      *info_filename = '\0';
      strcat (info_filename, outputprefix);
      strcpy (&info_filename[prefixlen], info_dens_label);
      strcpy (&info_filename[prefixlen + info_label_len], "\0");
      finfo = fopen (info_filename, "w");
      fprintf (finfo, "Grid order - row major order (C):\n x_0 y_0 z_0\n x_0 y_0 "
                      "z_1\n ...\n x_0 y_0 z_n-1\n x_0 y_1 "
                      "z_0\n ...\n\n");
      fprintf (finfo, "Index:   X coordinate:\n");
      for (i = 0; i < nx; i++)
        fprintf (finfo, " %4d       %9.4lf\n", i, xmin + (i + 0.5) * dx);
      fprintf (finfo, "\n");
      fprintf (finfo, "Index:   Y coordinate:\n");
      for (j = 0; j < ny; j++)
        fprintf (finfo, " %4d       %9.4lf\n", j, ymin + (j + 0.5) * dy);
      fprintf (finfo, "\n");
      fprintf (finfo, "Index:   Z coordinate:\n");
      for (k = 0; k < nz; k++)
        fprintf (finfo, " %4d       %9.4lf\n", k, zmin + (k + 0.5) * dz);
      fprintf (finfo, "\n");
      fprintf (finfo, "Contents of the output binary file: \n");
      fprintf (finfo, "Events (long int):    %ld\n", nevents);
      fprintf (finfo, "Time (double)\n");
      fprintf (finfo, "Output content flag (int):   %d\n", output_content_info);
      fprintf (finfo, "Hadrons species np (int):    %d\n", np);
      if (include_resonances)
        {
          fprintf (finfo, "Resonance species nr (int):    %d\n", nr);
        }
      fprintf (finfo, "x cells nx (int):    %d\n", nx);
      fprintf (finfo, "y cells ny (int):    %d\n", ny);
      fprintf (finfo, "z cells nz (int):    %d\n", nz);
      fprintf (finfo, "dx (double):    %lf\n", dx);
      fprintf (finfo, "dy (double):    %lf\n", dy);
      fprintf (finfo, "dz (double):    %lf\n", dz);
      fprintf (finfo, "xmin, i.e. low border of the 1st x cell (double):    %lf\n", xmin);
      fprintf (finfo, "ymin, i.e. low border of the 1st y cell (double):    %lf\n", ymin);
      fprintf (finfo, "zmin, i.e. low border of the 1st z cell (double):    %lf\n", zmin);
      fprintf (finfo, "output_content_info (int), optional output flag:    %d\n", output_content_info);
      fprintf (finfo,
               "Volume factor dxdydz (divide by this number to get a density): "
               "%lf\n\n",
               dx * dy * dz);
      fprintf (finfo, "Indexes for Tmunu (%d entries):\n", np);
      for (i = 0; i < np - 1; i++)
        {
          fprintf (finfo, "%4d:\t%s\n", i, associate_particle_array_index (i));
        }
      fprintf (finfo, "%4d:\tAll other particles\n", i);
      if (include_resonances)
        {
          fprintf (finfo, "Indexes for resonance Tmunu (%d entries):\n", nr);
          for (i = 0; i < nr; i++)
            {
              fprintf (finfo, "%4d:\t%s\n", i, associate_resonance_array_index (i));
            }
        }
      fclose (finfo);
      printf ("Written grid and particle informations in file %s\n", info_filename);
      free (info_filename);
      first_time = 1;
    }

  if (include_total_baryon)
    {
      offset += 1;
    }

  empty_array = (double *)calloc (offset, sizeof (double));
  if (empty_array == NULL)
    {
      printf ("Unable to create the empty_array in function write_densities in "
              "io.c. This is a weird situation and I "
              "prefer to quit.\n");
      exit (4);
    }

  for (h = 0; h < nt; h++)
    {
      if (h == 0)
        printf ("Now computing and printing densities of the total %ld events\n", nevents);
      snprintf (time_suffix, 8, "%07.3lf", time_int_array[h]);
      strncpy (&time_suffix[3], "_", 1); // we replace the dot with an underscore
      output_filename_Tp = (char *)malloc (sizeof (char) * (prefixlen + dens_label_len + 8));
      *output_filename_Tp = '\0';
      strncat (output_filename_Tp, outputprefix, prefixlen);
      strncat (output_filename_Tp, dens_label, dens_label_len);
      strncat (output_filename_Tp, time_suffix, 8);

      fTp = fopen (output_filename_Tp, "wb");
      if (fTp == NULL)
        {
          printf ("Sorry, I could not open the output file %s, so I am forced to "
                  "quit.\n",
                  output_filename_Tp);
          exit (3);
        }

      fwrite (&nevents, sizeof (long int), 1, fTp);
#ifdef DBG_DENS
      printf ("nevents: %ld\n", nevents);
#endif
      fwrite (&time_int_array[h], sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("time: %lf\n", time_int_array[h]);
#endif
      fwrite (&output_content_info, sizeof (int), 1, fTp);
#ifdef DBG_DENS
      printf ("content info: %d\n", output_content_info);
#endif
      fwrite (&np, sizeof (int), 1, fTp);
#ifdef DBG_DENS
      printf ("np: %d\n", np);
#endif
      if (include_resonances)
        {
          fwrite (&nr, sizeof (int), 1, fTp);
#ifdef DBG_DENS
          printf ("nr: %d\n", nr);
#endif
        }
      fwrite (&nx, sizeof (int), 1, fTp);
#ifdef DBG_DENS
      printf ("nx: %d\n", nx);
#endif
      fwrite (&ny, sizeof (int), 1, fTp);
#ifdef DBG_DENS
      printf ("ny: %d\n", ny);
#endif
      fwrite (&nz, sizeof (int), 1, fTp);
#ifdef DBG_DENS
      printf ("nz: %d\n", nz);
#endif
      fwrite (&dx, sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("dx: %lf\n", dx);
#endif
      fwrite (&dy, sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("dy: %lf\n", dy);
#endif
      fwrite (&dz, sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("dz: %lf\n", dz);
#endif
      fwrite (&xmin, sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("xmin: %lf\n", xmin);
#endif
      fwrite (&ymin, sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("ymin: %lf\n", ymin);
#endif
      fwrite (&zmin, sizeof (double), 1, fTp);
#ifdef DBG_DENS
      printf ("zmin: %lf\n", zmin);
      int w = 0;
#endif

      for (i = 0; i < nx; i++)
        {
          for (j = 0; j < ny; j++)
            {
              for (k = 0; k < nz; k++)
                {

                  jf_arg = (Jb[J0 + JBL] * Jb[J0 + JBL] - Jb[J1 + JBL] * Jb[J1 + JBL] - Jb[J2 + JBL] * Jb[J2 + JBL]
                            - Jb[J3 + JBL] * Jb[J3 + JBL]);
                  if (jf_arg > 0)
                    {
                      jf = sqrt (jf_arg);
                      u4[0] = Jb[J0 + JBL] / jf;
                      u4[1] = Jb[J1 + JBL] / jf;
                      u4[2] = Jb[J2 + JBL] / jf;
                      u4[3] = Jb[J3 + JBL] / jf;
                      // printf("*** Point with coordinates: %d  %d  %d  %e  %e
                      // %e\n",i,j,k,xmin+(i+0.5)*dx,ymin+(j+0.5)*dy,zmin+(k+0.5)*dz);
                      // printf("*** Point with coordinates: %e  %e
                      // %e\n",xmin+(i+0.5)*dx,ymin+(j+0.5)*dy,zmin+(k+0.5)*dz); printf("u
                      // components: %e  %e  %e  %e\n",u4[0], u4[1], u4[2], u4[3]);
                      rho = (Jb[J0 + JBL] * u4[0] - Jb[J1 + JBL] * u4[1] - Jb[J2 + JBL] * u4[2] - Jb[J3 + JBL] * u4[3])
                            / (cell_volume * nevents);
                      // printf("Jb components and rho: %e  %e  %e  %e  %e\n",Jb[J0+JBL],
                      // Jb[J1+JBL], Jb[J2+JBL], Jb[J3+JBL], rho);
                      fwrite (&rho, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                      printf ("w, i, j, k, rho: %d  %d  %d  %d  %lf\n", w, i, j, k, rho);
                      w++;
#endif
                      for (l = 1; l < 4; l++)
                        vel_B[l - 1] = u4[l] / u4[0];
                      fwrite (vel_B, sizeof (double), 3, fTp);
#ifdef DBG_DENS
                      printf ("w, i, j, k, vel_B: %d  %d  %d  %d  %lf  %lf   %lf\n", w, i, j, k, vel_B[0], vel_B[1],
                              vel_B[2]);
                      w += 3;
#endif
                      // printf("vel components: %lf  %lf  %lf  \n",vel_B[0],vel_B[1],vel_B[2]);
                      if (include_total_baryon)
                        {
                          rho_t = (Jt[J0 + JBL] * u4[0] - Jt[J1 + JBL] * u4[1] - Jt[J2 + JBL] * u4[2]
                                   - Jt[J3 + JBL] * u4[3])
                                  / (cell_volume * nevents);
                          fwrite (&rho_t, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                          printf ("w, i, j, k, rho tot bar: %d  %d  %d  %d  %lf\n", w, i, j, k, rho_t);
                          w++;
#endif
                        }
                      rho_c
                          = (Jc[J0 + JBL] * u4[0] - Jc[J1 + JBL] * u4[1] - Jc[J2 + JBL] * u4[2] - Jc[J3 + JBL] * u4[3])
                            / (cell_volume * nevents);
                      rho_s
                          = (Js[J0 + JBL] * u4[0] - Js[J1 + JBL] * u4[1] - Js[J2 + JBL] * u4[2] - Js[J3 + JBL] * u4[3])
                            / (cell_volume * nevents);
                      fwrite (&rho_c, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                      printf ("w, i, j, k, rho c: %d  %d  %d  %d  %lf\n", w, i, j, k, rho_c);
                      w++;
#endif
                      fwrite (&rho_s, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                      printf ("w,i, j, k, rho s: %d  %d  %d  %d  %lf\n", w, i, j, k, rho_s);
                      w++;
#endif
                      // printf("Jt components and rho_t: %e  %e  %e  %e %e\n",Jt[J0+JBL],
                      // Jt[J1+JBL], Jt[J2+JBL], Jt[J3+JBL], rho_t); printf("Jc components
                      // and rho_c: %e  %e  %e  %e  %e\n",Jc[J0+JBL], Jc[J1+JBL],
                      // Jc[J2+JBL], Jc[J3+JBL], rho_c); printf("Js components and rho_s:
                      // %e  %e  %e  %e %e\n",Js[J0+JBL], Js[J1+JBL], Js[J2+JBL],
                      // Js[J3+JBL], rho_s);
                      Ic_diffcheck[0] = (1 - u4[0] * u4[0]) * Jc[J0 + JBL] + u4[0] * u4[1] * Jc[J1 + JBL]
                                        + u4[0] * u4[2] * Jc[J2 + JBL] + u4[0] * u4[3] * Jc[J3 + JBL];
                      Ic_diffcheck[1] = -u4[0] * u4[1] * Jc[J0 + JBL] + (1 + u4[1] * u4[1]) * Jc[J1 + JBL]
                                        + u4[1] * u4[2] * Jc[J2 + JBL] + u4[1] * u4[3] * Jc[J3 + JBL];
                      Ic_diffcheck[2] = -u4[0] * u4[2] * Jc[J0 + JBL] + u4[2] * u4[1] * Jc[J1 + JBL]
                                        + (1 + u4[2] * u4[2]) * Jc[J2 + JBL] + u4[2] * u4[3] * Jc[J3 + JBL];
                      Ic_diffcheck[3] = -u4[0] * u4[3] * Jc[J0 + JBL] + u4[3] * u4[1] * Jc[J1 + JBL]
                                        + u4[3] * u4[2] * Jc[J2 + JBL] + (1 + u4[3] * u4[3]) * Jc[J3 + JBL];
                      Is_diffcheck[0] = (1 - u4[0] * u4[0]) * Js[J0 + JBL] + u4[0] * u4[1] * Js[J1 + JBL]
                                        + u4[0] * u4[2] * Js[J2 + JBL] + u4[0] * u4[3] * Js[J3 + JBL];
                      Is_diffcheck[1] = -u4[0] * u4[1] * Js[J0 + JBL] + (1 + u4[1] * u4[1]) * Js[J1 + JBL]
                                        + u4[1] * u4[2] * Js[J2 + JBL] + u4[1] * u4[3] * Js[J3 + JBL];
                      Is_diffcheck[2] = -u4[0] * u4[2] * Js[J0 + JBL] + u4[2] * u4[1] * Js[J1 + JBL]
                                        + (1 + u4[2] * u4[2]) * Js[J2 + JBL] + u4[2] * u4[3] * Js[J3 + JBL];
                      Is_diffcheck[3] = -u4[0] * u4[3] * Js[J0 + JBL] + u4[3] * u4[1] * Js[J1 + JBL]
                                        + u4[3] * u4[2] * Js[J2 + JBL] + (1 + u4[3] * u4[3]) * Js[J3 + JBL];
                      for (l = 0; l < 4; l++)
                        Ic_diffusion[l] = Jc[l + JBL] / (cell_volume * nevents) - rho_c * u4[l];
                      for (l = 0; l < 4; l++)
                        Is_diffusion[l] = Js[l + JBL] / (cell_volume * nevents) - rho_s * u4[l];
                      for (l = 0; l < 4; l++)
                        {
                          if (fabs (Ic_diffusion[l] - Ic_diffcheck[l] / (cell_volume * nevents)) > 1.e-10)
                            printf ("Warning, mismatching in charge diffusion currents at "
                                    "i=%d, j=%d, k=%d! %lf vs %lf\n",
                                    i, j, k, Ic_diffusion[l], Ic_diffcheck[l]);
                          if (fabs (Is_diffusion[l] - Is_diffcheck[l] / (cell_volume * nevents)) > 1.e-10)
                            printf ("Warning, mismatching in strange diffusion currents at "
                                    "i=%d, j=%d, k=%d! %lf vs %lf\n",
                                    i, j, k, Is_diffusion[l], Is_diffcheck[l]);
                        }
                      fwrite (Ic_diffusion, sizeof (double), 4, fTp);
#ifdef DBG_DENS
                      printf ("w, i, j, k, Ic components: %d  %d  %d  %d  %lf  %lf  %lf  %lf\n", w, i, j, k,
                              Ic_diffusion[0], Ic_diffusion[1], Ic_diffusion[2], Ic_diffusion[3]);
                      w += 4;
#endif
                      fwrite (Is_diffusion, sizeof (double), 4, fTp);
#ifdef DBG_DENS
                      printf ("w, i, j, k, Is components: %d  %d  %d  %d  %lf  %lf  %lf  %lf\n", w, i, j, k,
                              Is_diffusion[0], Is_diffusion[1], Is_diffusion[2], Is_diffusion[3]);
                      w += 4;
#endif
                      // printf("Ic diffusion components: %e  %e  %e
                      // %e\n",Ic_diffusion[0],Ic_diffusion[1],Ic_diffusion[2],
                      // Ic_diffusion[3]); printf("Is diffusion components: %e  %e  %e
                      // %e\n",Is_diffusion[0],Is_diffusion[1],Is_diffusion[2],
                      // Is_diffusion[3]); printf("\n\n");
                    }
                  else
                    {
                      if (jf_arg < 0)
                        {
                          printf ("Problem with Jb at index h: %d, i: %d, j: %d, k: %d\n", h, i, j, k);
                          printf ("Jb components: %9.5e    %9.5e    %9.5e    %9.5e\n", Jb[J0 + JBL], Jb[J1 + JBL],
                                  Jb[J2 + JBL], Jb[J3 + JBL]);
                        }
                      fwrite (empty_array, sizeof (double), offset,
                              fTp); // it includes rho,velB,rho_t,rho_c, rho_s,
                                    // Ic_diffusion,Is_diffusion, rho and eps for all hadrons
                    }
                  for (p = 0; p < np; p++)
                    {
                      if (jf_arg > 0)
                        {
                          rho = (Jp[J0 + JPL] * u4[0] - Jp[J1 + JPL] * u4[1] - Jp[J2 + JPL] * u4[2]
                                 - Jp[J3 + JPL] * u4[3])
                                / (cell_volume * nevents);
                          eps = ((u4[0] * Tp[T00 + TLOC] - u4[1] * Tp[T10 + TLOC] - u4[2] * Tp[T20 + TLOC]
                                  - u4[3] * Tp[T30 + TLOC])
                                     * u4[0]
                                 - (u4[0] * Tp[T01 + TLOC] - u4[1] * Tp[T11 + TLOC] - u4[2] * Tp[T21 + TLOC]
                                    - u4[3] * Tp[T31 + TLOC])
                                       * u4[1]
                                 - (u4[0] * Tp[T02 + TLOC] - u4[1] * Tp[T12 + TLOC] - u4[2] * Tp[T22 + TLOC]
                                    - u4[3] * Tp[T32 + TLOC])
                                       * u4[2]
                                 - (u4[0] * Tp[T03 + TLOC] - u4[1] * Tp[T13 + TLOC] - u4[2] * Tp[T23 + TLOC]
                                    - u4[3] * Tp[T33 + TLOC])
                                       * u4[3])
                                / (cell_volume * nevents);
                        }
                      jh_arg = (Jp[J0 + JPL] * Jp[J0 + JPL] - Jp[J1 + JPL] * Jp[J1 + JPL] - Jp[J2 + JPL] * Jp[J2 + JPL]
                                - Jp[J3 + JPL] * Jp[J3 + JPL]);
                      if (jh_arg > 0)
                        {
                          jf = sqrt (jh_arg);
                          u4[0] = Jp[J0 + JPL] / jf;
                          u4[1] = Jp[J1 + JPL] / jf;
                          u4[2] = Jp[J2 + JPL] / jf;
                          u4[3] = Jp[J3 + JPL] / jf;
                          for (l = 1; l < 4; l++)
                            vel_h[l - 1] = u4[l] / u4[0];
                        }
                      tmp_value = (double)Pnum[PNLOC];
                      fwrite (&tmp_value, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                      printf ("w, i, j, k, p, n p: %d  %d  %d  %d  %d  %lf\n", w, i, j, k, p, tmp_value);
                      w++;
#endif
                      if (jf_arg > 0)
                        {
                          fwrite (&rho, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                          printf ("w, i, j, k, p, rho p: %d  %d  %d  %d  %d  %lf\n", w, i, j, k, p, rho);
                          w++;
#endif
                          fwrite (&eps, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                          printf ("w, i, j, k, p, eps p: %d  %d  %d  %d  %d  %lf\n", w, i, j, k, p, eps);
                          w++;
#endif
                        }
                      else
                        {
                          fwrite (&null_2, sizeof (double), 2, fTp);
#ifdef DBG_DENS
                          w += 2;
#endif
                        }
                      if (jh_arg > 0)
                        {
                          fwrite (vel_h, sizeof (double), 3, fTp);
#ifdef DBG_DENS
                          printf ("w, i, j, k, p, vel p: %d  %d  %d  %d  %d  %lf  %lf  %lf\n", w, i, j, k, p, vel_h[0],
                                  vel_h[1], vel_h[2]);
#endif
                        }
                      else
                        {
                          fwrite (vel_null, sizeof (double), 3, fTp);
                        }
#ifdef DBG_DENS
                      w += 3;
#endif
                    } // end cycle over p
                  if (include_resonances)
                    {
                      for (r = 0; r < nr; r++)
                        {
                          if (jf_arg > 0)
                            {
                              rho = (Jr[J0 + JRL] * u4[0] - Jr[J1 + JRL] * u4[1] - Jr[J2 + JRL] * u4[2]
                                     - Jr[J3 + JRL] * u4[3])
                                    / (cell_volume * nevents);
                              eps = ((u4[0] * Tr[T00 + RTLOC] - u4[1] * Tr[T10 + RTLOC] - u4[2] * Tr[T20 + RTLOC]
                                      - u4[3] * Tr[T30 + RTLOC])
                                         * u4[0]
                                     - (u4[0] * Tr[T01 + RTLOC] - u4[1] * Tr[T11 + RTLOC] - u4[2] * Tr[T21 + RTLOC]
                                        - u4[3] * Tr[T31 + RTLOC])
                                           * u4[1]
                                     - (u4[0] * Tr[T02 + RTLOC] - u4[1] * Tr[T12 + RTLOC] - u4[2] * Tr[T22 + RTLOC]
                                        - u4[3] * Tr[T32 + RTLOC])
                                           * u4[2]
                                     - (u4[0] * Tr[T03 + RTLOC] - u4[1] * Tr[T13 + RTLOC] - u4[2] * Tr[T23 + RTLOC]
                                        - u4[3] * Tr[T33 + RTLOC])
                                           * u4[3])
                                    / (cell_volume * nevents);
                            }
                          jh_arg = (Jr[J0 + JRL] * Jr[J0 + JRL] - Jr[J1 + JRL] * Jr[J1 + JRL]
                                    - Jr[J2 + JRL] * Jr[J2 + JRL] - Jr[J3 + JRL] * Jr[J3 + JRL]);
                          if (jh_arg > 0)
                            {
                              jf = sqrt (jh_arg);
                              u4[0] = Jr[J0 + JRL] / jf;
                              u4[1] = Jr[J1 + JRL] / jf;
                              u4[2] = Jr[J2 + JRL] / jf;
                              u4[3] = Jr[J3 + JRL] / jf;
                              for (l = 1; l < 4; l++)
                                vel_h[l - 1] = u4[l] / u4[0];
                            }
                          tmp_value = (double)Rnum[RNLOC];
                          fwrite (&tmp_value, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                          printf ("w, i, j, k, r, n p: %d  %d  %d  %d  %d  %lf\n", w, i, j, k, r, tmp_value);
                          w++;
#endif
                          if (jf_arg > 0)
                            {
                              fwrite (&rho, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                              printf ("w, i, j, k, r, rho r: %d  %d  %d  %d  %d  %lf\n", w, i, j, k, r, rho);
                              w++;
#endif
                              fwrite (&eps, sizeof (double), 1, fTp);
#ifdef DBG_DENS
                              printf ("w, i, j, k, r, eps r: %d  %d  %d  %d  %d  %lf\n", w, i, j, k, r, eps);
                              w++;
#endif
                            }
                          else
                            {
                              fwrite (&null_2, sizeof (double), 2, fTp);
#ifdef DBG_DENS
                              w += 2;
#endif
                            }
                          if (jh_arg > 0)
                            {
                              fwrite (vel_h, sizeof (double), 3, fTp);
#ifdef DBG_DENS
                              printf ("w, i, j, k, p, vel r: %d  %d  %d  %d  %d  %lf  %lf  %lf\n", w, i, j, k, r,
                                      vel_h[0], vel_h[1], vel_h[2]);
#endif
                            }
                          else
                            {
                              fwrite (vel_null, sizeof (double), 3, fTp);
                            }
#ifdef DBG_DENS
                          w += 3;
#endif
                        }
                    }
                }
            }
        }

      printf ("Densities data saved in file %s.\n", output_filename_Tp);
      fflush (fTp);
      fclose (fTp);
      free (output_filename_Tp);
    }

  free (empty_array);
}

void
read_data (char *inputfile, double *data_Tp, double *data_Jp, double *data_Jb, double *data_Jc, double *data_Js,
           double *data_Jt, double *data_Jr, double *data_Tr, long int *data_Pnum, long int *data_Rnum,
           long int *nevents, double *data_time)
{
  FILE *infile;
  int np2, nr2, nx2, ny2, nz2;
  int output_content_info2;
  double dx2, dy2, dz2;
  double xmin2, ymin2, zmin2;

  infile = fopen (inputfile, "rb");
  if (infile == NULL)
    {
      printf ("Sorry, I could not open the input file %s, so I am forced to quit.\n", inputfile);
      exit (3);
    }
  ret_it = fread (nevents, sizeof (long int), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (data_time, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&output_content_info2, sizeof (int), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&np2, sizeof (int), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&nx2, sizeof (int), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&ny2, sizeof (int), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&nz2, sizeof (int), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dx2, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dy2, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dz2, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&xmin2, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&ymin2, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&zmin2, sizeof (double), 1, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  // at some point these checks should be changed and the comparison should be
  // only with the values in the first file, not with those hardcoded
  if (output_content_info != output_content_info2)
    {
      printf ("Error: mismatching between the output_content_info flags!!!\n");
      printf ("The hardcoded value is: output_content_info=%d\n", output_content_info);
      printf ("The new value in %s is: output_content_info2=%d\n", inputfile, output_content_info2);
      printf ("Sorry, but I have to quit..\n");
      fclose (infile);
      exit (6);
    }
  if (np != np2)
    {
      printf ("Error: mismatching between the number of particles!!!\n");
      printf ("The hardcoded value is: np=%d\n", np);
      printf ("The new value in %s is: np2=%d\n", inputfile, np2);
      printf ("Sorry, but I have to quit..\n");
      fclose (infile);
      exit (6);
    }
  if ((nx != nx2) || (ny != ny2) || (nz != nz2) || (dx != dx2) || (dy != dy2) || (dz != dz2) || (xmin != xmin2)
      || (ymin != ymin2) || (zmin != zmin2))
    {
      printf ("Error: mismatching between the dimension and/or resolution of the "
              "old and the new grid!!!\n");
      printf ("The old values are: nx=%d, ny=%d, nz=%d, dx=%6.4lf, dy=%6.4lf, "
              "dz=%6.4lf, xmin=%6.4lf, ymin=%6.4lf, "
              "zmin=%6.4lf\n",
              nx, ny, nz, dx, dy, dz, xmin, ymin, zmin);
      printf ("The new values in %s are: nx2=%d, ny2=%d, nz2=%d, dx2=%6.4lf, "
              "dy2=%6.4lf, dz2=%6.4lf, xmin2=%6.4lf, "
              "ymin2=%6.4lf, zmin2=%6.4lf\n",
              inputfile, nx2, ny2, nz2, dx2, dy2, dz2, xmin2, ymin2, zmin2);
      printf ("Sorry, but I have to quit..\n");
      fclose (infile);
      exit (6);
    }
  ret_it = fread (data_Tp, sizeof (double), nx * ny * nz * np * 10, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (data_Jp, sizeof (double), nx * ny * nz * np * 4, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (data_Jb, sizeof (double), nx * ny * nz * 4, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (data_Jc, sizeof (double), nx * ny * nz * 4, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (data_Js, sizeof (double), nx * ny * nz * 4, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  if (include_total_baryon)
    {
      ret_it = fread (data_Jt, sizeof (double), nx * ny * nz * 4, infile);
      if (ret_it == 0)
        {
          printf ("Failure in reading data. Exiting.\n");
          exit (4);
        }
    }
  if (include_resonances)
    {
      ret_it = fread (&nr2, sizeof (int), 1, infile);
      if (ret_it == 0)
        {
          printf ("Failure in reading data. Exiting.\n");
          exit (4);
        }
      else
        {
          if (nr != nr2)
            {
              printf ("Error: mismatching between the number of resonances!!!\n");
              printf ("The hardcoded value is: nr=%d\n", nr);
              printf ("The new value in %s is: nr2=%d\n", inputfile, nr2);
              printf ("Sorry, but I have to quit..\n");
              fclose (infile);
              exit (6);
            }
        }
      ret_it = fread (data_Jr, sizeof (double), nx * ny * nz * nr * 4, infile);
      if (ret_it == 0)
        {
          printf ("Failure in reading data. Exiting.\n");
          exit (4);
        }
      ret_it = fread (data_Tr, sizeof (double), nx * ny * nz * nr * 10, infile);
      if (ret_it == 0)
        {
          printf ("Failure in reading data. Exiting.\n");
          exit (4);
        }
    }
  ret_it = fread (data_Pnum, sizeof (long int), nx * ny * nz * np, infile);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  if (include_resonances)
    {
      ret_it = fread (data_Rnum, sizeof (long int), nx * ny * nz * nr, infile);
      if (ret_it == 0)
        {
          printf ("Failure in reading data. Exiting.\n");
          exit (4);
        }
    }
  printf ("%s read.\n", inputfile);
  fclose (infile);
}
