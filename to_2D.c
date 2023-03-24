#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
help ()
{
  printf (
      "Syntax: ./to_2D.exe <inputfile> <outputfile> <i,j or k> <index>\nIt prints a 2D slice of the energy density,"
      " suitable to be plotted with gnuplot, at fixed i, j or k equal to index.\nThe informations about the grid"
      " and the particles are contained into definitions.h.\n\n");
}

int
main (int argc, char *argv[])
{

  FILE *fin, *fout;
  long int nevents;
  double time;
  int i, j, k, p, l, chosen_index, chosen_dir, inner_counter;
  double nump, bb, had_edens, edens, xmin, ymin, zmin, density;
  double bb3[3];
  double *bitbucket;
  int nx, ny, nz, np;
  double dx, dy, dz;
  size_t ret_it;
  const int entries = 6;

  if (argc != 5)
    {
      help ();
      exit (1);
    }

  if (strncmp (argv[3], "i", 1) == 0)
    {
      chosen_dir = 0;
    }
  else if (strncmp (argv[3], "j", 1) == 0)
    {
      chosen_dir = 1;
    }
  else if (strncmp (argv[3], "k", 1) == 0)
    {
      chosen_dir = 2;
    }
  else
    {
      printf ("Sorry, unable to determine the orientation of the slicing plane\n");
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
      fclose (fin);
      exit (2);
    }

  chosen_index = atoi (argv[4]);

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

  bitbucket = (double *)malloc (sizeof (double) * (15 + entries * np));
  if (bitbucket == NULL)
    {
      printf ("Allocation of bitbucket array failed. I quit.\n");
      exit (2);
    }
  inner_counter = 0;
  for (i = 0; i < nx; i++)
    {
      for (j = 0; j < ny; j++)
        {
          for (k = 0; k < nz; k++)
            {
              if (((chosen_dir == 0) && (i == chosen_index)) || ((chosen_dir == 1) && (j == chosen_index))
                  || ((chosen_dir == 2) && (k == chosen_index)))
                {
                  fprintf (fout, "%7.3e  %7.3e  %7.3e  ", xmin + (i + 0.5) * dx, ymin + (j + 0.5) * dy,
                           zmin + (k + 0.5) * dz);
                  for (l = 0; l < 15; l++)
                    {
                      ret_it = fread (&density, sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      fprintf (fout, "%10.6e  ", density);
                    }
                  edens = 0;
                  for (p = 0; p < np; p++)
                    {
                      ret_it = fread (&nump, sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      ret_it = fread (&bb, sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      ret_it = fread (&had_edens, sizeof (double), 1, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      ret_it = fread (bb3, sizeof (double), 3, fin);
                      if (ret_it == 0)
                        {
                          printf ("Failure in reading data. Exiting.\n");
                          exit (4);
                        }
                      edens = edens + had_edens;
                    }
                  fprintf (fout, "%10.6e  ", edens); // total energy density
                  fprintf (fout, "\n");
                  inner_counter = inner_counter + 1;
                  // this counter tells when to add an additional new line at the end of a block
                  // please, note that for slices along i and j it is nz, while it is ny for slices along k
                  if (((chosen_dir == 0) && (inner_counter == nz)) || ((chosen_dir == 1) && (inner_counter == nz))
                      || ((chosen_dir == 2) && (inner_counter == ny)))
                    {
                      inner_counter = 0;
                      fprintf (fout, "\n");
                    }
                }
              else
                {
                  ret_it = fread (bitbucket, sizeof (double), 15 + entries * np, fin);
                  if (ret_it == 0)
                    {
                      printf ("Failure in reading data. Exiting.\n");
                      exit (4);
                    }
                }
            }
        }
    }

  free (bitbucket);
  fclose (fin);
  fclose (fout);
  return 0;
}
