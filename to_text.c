#include <stdio.h>
#include <stdlib.h>

void
help ()
{
  printf ("Syntax: ./to_text.exe <inputfile> <outputfile>\nThe informations about the grid and the particles are "
          "contained into definitions.h.\n\n");
}

int
main (int argc, char *argv[])
{

  FILE *fin, *fout;
  double *datap;
  long int nevents;
  double time, xmin, ymin, zmin;
  int i, j, k, p, h, l;
  int nx, ny, nz, np;
  double dx, dy, dz;
  double entot;
  size_t ret_it;
  const int p_entries = 6;

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
  ret_it = fread (&np, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data. Exiting.\n");
      exit (4);
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

  datap = (double *)malloc (sizeof (double) * ny * nz * (p_entries * np + 15));

  if (datap == NULL)
    {
      printf ("Sorry, but I cannot allocate the datap array to temporary store the input data...\n");
      exit (3);
    }

  fprintf (fout, "Number of events: %12ld  \n", nevents);
  fprintf (fout, "Time: %12.3e  \n", time);
  fprintf (fout, "Number of particles: %12d  \n", np);
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
              ret_it = fread (&datap[h], sizeof (double), 15, fin); // rho,vx,vy,vz,rho_c,rho_s,rho_t,Ic_diff,Is_diff
              if (ret_it == 0)
                {
                  printf ("Failure in reading data. Exiting.\n");
                  exit (4);
                }
              if (ret_it == 0)
                {
                  printf ("Failure in reading data. Exiting.\n");
                  exit (4);
                }
              h += 15;
              for (p = 0; p < np; p++)
                {
                  ret_it = fread (&datap[h + p_entries * p], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + p_entries * p + 1], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + p_entries * p + 2], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + p_entries * p + 3], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + p_entries * p + 4], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                  ret_it = fread (&datap[h + p_entries * p + 5], sizeof (double), 1, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                }
              h += p_entries * p;
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
              fprintf (fout, "total baryon + antibaryon density: %14.9e\n", datap[h + 4]);
              fprintf (fout, "charge density: %14.9e\n", datap[h + 5]);
              fprintf (fout, "strangeness density: %14.9e\n", datap[h + 6]);
              fprintf (fout, "electric charge diffusion current components: %14.9e, %14.9e, %14.9e, %14.9e\n",
                       datap[h + 7], datap[h + 8], datap[h + 9], datap[h + 10]);
              fprintf (fout, "strangeness diffusion current components: %14.9e, %14.9e, %14.9e, %14.9e\n",
                       datap[h + 11], datap[h + 12], datap[h + 13], datap[h + 14]);
              h += 15;
              entot = 0.;
              for (p = 0; p < np; p++)
                {
                  entot = entot + datap[p_entries * p + h + 2];
                  fprintf (
                      fout,
                      "Hadron kind: %4d, total cell number: %14ld, number density: %14.9e, energy density: %14.9e,\
                       vx: %14.9e, vy: %14.9e, vz: %14.9e\n",
                      p, (long int)datap[p_entries * p + h], datap[p_entries * p + h + 1],
                      datap[p_entries * p + h + 2], datap[p_entries * p + h + 3], datap[p_entries * p + h + 4],
                      datap[p_entries * p + h + 5]);
                }
              fprintf (fout, "total energy density: %14.9e\n", entot);
              fprintf (fout, "\n");
              h += p_entries * p;
            }
        }
    }

  free (datap);
  fclose (fout);
  fclose (fin);
  return 0;
}
