/* Author: Gabriele Inghirami g.inghirami@gsi.de (2019-2022) - License: GPLv3 */

#include <stdio.h>
#include <stdlib.h>

const int shift_total_baryon_on = 10; // value added to output_content_info if
                                      // the total baryons are considered
const int shift_resonances_on = 100;  // value added to output_content_info if the resonances are considered
int total_baryon_included = 0;        // flag that will be set to > 0 if total_baryons
                                      // are included, to < 0 if not
int resonances_included = 0;          // flag that will be set to > 0 if resonances are included, to < 0 if not
size_t ret_it;

void
help ()
{
  printf ("Syntax: ./to_text.exe <inputfile> <outputfile>\nThe informations "
          "about the grid and the particles are "
          "contained into definitions.h.\n\n");
}

void
interpret_output_content (int num)
{
  if (num < shift_resonances_on)
    {
      resonances_included = -1;
    }
  else
    {
      resonances_included = 1;
      num -= resonances_included;
    }
  if (num < shift_total_baryon_on)
    {
      total_baryon_included = -1;
    }
  else
    {
      total_baryon_included = 1;
    }
}

int
main (int argc, char *argv[])
{

  FILE *fin, *fout;
  double *datap;
  long int nevents;
  double time, xmin, ymin, zmin;
  int i, j, k, p, h, l, r;
  int nx, ny, nz, np, nr;
  double dx, dy, dz;
  double entot;
  int output_content;
  int offset = 14;
  int offset_res = 0;

  if (argc != 3)
    {
      help ();
      exit (1);
    }

  fin = fopen (argv[1], "r");
  if (fin == NULL)
    {
      printf ("Sorry, I was unable to open the input file %s\n", argv[1]);
      exit (2);
    }

  fout = fopen (argv[2], "w+");
  if (fout == NULL)
    {
      printf ("Sorry, I was unable to create the output file %s\n", argv[2]);
      exit (2);
    }
  ret_it = fread (&nevents, sizeof (long int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&time, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&output_content, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  interpret_output_content (output_content);
  ret_it = fread (&np, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  if (resonances_included > 0)
    {
      ret_it = fread (&nr, sizeof (int), 1, fin);
      if (ret_it == 0)
        {
          printf ("Failure in reading data. Exiting.\n");
          exit (4);
        }
    }
  ret_it = fread (&nx, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&ny, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&nz, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dx, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dy, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dz, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&xmin, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&ymin, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }
  ret_it = fread (&zmin, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
    }

  if (total_baryon_included > 0)
    offset += 1;
  if (resonances_included > 0)
    offset_res += 3 * nr;

  datap = (double *)malloc (sizeof (double) * ny * nz * (3 * np + offset + offset_res));

  if (datap == NULL)
    {
      printf ("Sorry, but I cannot allocate the datap array to temporary store "
              "the input data...\n");
      exit (3);
    }

  fprintf (fout, "Number of events: %12ld  \n", nevents);
  fprintf (fout, "Time: %12.3e  \n", time);
  fprintf (fout, "Number of particles: %12d  \n", np);
  if (resonances_included > 0)
    fprintf (fout, "Number of resonances: %12d  \n", nr);
  fprintf (fout, "nx: %6d, ny: %6d, nz: %6d\n", nx, ny, nz);
  fprintf (fout, "dx: %7.3e, dy: %7.3e, dz: %7.3e\n", dx, dy, dz);
  fprintf (fout, "xmin: %7.3e, ymin: %7.3e, zmin: %7.3e\n", dx, dy, dz);

  for (i = 0; i < nx; i++)
    {
      h = 0;
      for (j = 0; j < ny; j++)
        {
          for (k = 0; k < nz; k++)
            {
              ret_it = fread (&datap[h], sizeof (double), offset,
                              fin); // rho,vx,vy,vz,rho_c,rho_s,rho_t,Ic_diff,Is_diff
              if (ret_it == 0)
                {
                  printf ("Failure in reading data. Exiting.\n");
                  exit (4);
                }
              h += offset;
              for (p = 0; p < np; p++)
                {
                  ret_it = fread (&datap[h + 3 * p], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + 3 * p + 1], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + 3 * p + 2], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                }
              h += 3 * p;
              if (resonances_included > 0)
                {
                  for (r = 0; r < nr; r++)
                    {
                      ret_it = fread (&datap[h + 3 * r], sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      ret_it = fread (&datap[h + 3 * r + 1], sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      ret_it = fread (&datap[h + 3 * r + 2], sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                    }
                  h += 3 * r;
                }
            }
        }
      h = 0;
      for (j = 0; j < ny; j++)
        {
          for (k = 0; k < nz; k++)
            {
              fprintf (fout, "Position: %7.3e  %7.3e  %7.3e  \n", xmin + (i + 0.5) * dx, ymin + (j + 0.5) * dy,
                       zmin + (k + 0.5) * dz);
              fprintf (fout, "baryon density: %14.9e\n", datap[h]);
              fprintf (fout, "vx: %14.9e, vy: %14.9e, vz: %14.9e\n", datap[h + 1], datap[h + 2], datap[h + 3]);
              h += 4;
              if (total_baryon_included > 0)
                {
                  fprintf (fout, "total baryon + antibaryon density: %14.9e\n", datap[h]);
                  h += 1;
                }
              fprintf (fout, "charge density: %14.9e\n", datap[h]);
              h += 1;
              fprintf (fout, "strangeness density: %14.9e\n", datap[h]);
              h += 1;
              fprintf (fout,
                       "electric charge diffusion current components: %14.9e, %14.9e, "
                       "%14.9e, %14.9e\n",
                       datap[h], datap[h + 1], datap[h + 2], datap[h + 3]);
              h += 4;
              fprintf (fout,
                       "strangeness diffusion current components: %14.9e, %14.9e, "
                       "%14.9e, %14.9e\n",
                       datap[h], datap[h + 1], datap[h + 2], datap[h + 3]);
              h += 4;
              entot = 0.;
              for (p = 0; p < np; p++)
                {
                  entot = entot + datap[3 * p + h + 2];
                  fprintf (fout,
                           "Hadron kind: %4d, total cell number: %14ld, number density: "
                           "%14.9e, energy density: %14.9e\n",
                           p, (long int)datap[3 * p + h], datap[3 * p + h + 1], datap[3 * p + h + 2]);
                }
              fprintf (fout, "total energy density: %14.9e\n", entot);
              fprintf (fout, "\n");
              h += 3 * p;
              if (resonances_included > 0)
                {
                  for (r = 0; r < nr; r++)
                    {
                      fprintf (fout,
                               "Resonance kind: %4d, total cell number: %14ld, number "
                               "density: %14.9e, energy "
                               "density: %14.9e\n",
                               r, (long int)datap[3 * r + h], datap[3 * r + h + 1], datap[3 * r + h + 2]);
                    }
                  fprintf (fout, "\n");
                  h += 3 * r;
                }
            }
        }
    }

  free (datap);
  fclose (fout);
  fclose (fin);
  return 0;
}
