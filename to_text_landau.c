#include <stdio.h>
#include <stdlib.h>

void
help ()
{
  printf ("Syntax: ./to_text_landau_density.exe <density_inputfile> <density_outputfile>\n\n");
}

int
main (int argc, char *argv[])
{

  FILE *fin, *fout;
  double *data;
  long int nevents;
  double time, xmin, ymin, zmin;
  int i, j, k, p, h, l;
  int nx, ny, nz, np;
  int bp;
  double dx, dy, dz;
  double Pnum, rho, eps;
  size_t ret_it;
  double u4[4], Ib_diffusion, Is_diffusion, Ic_diffusion;
  double TmunuL[10];
  

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
      printf ("Failure in reading data (nevents). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&time, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (time). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&np, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (np). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&nx, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (nx). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&ny, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (ny). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&nz, sizeof (int), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (nz). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dx, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (dx). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dy, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (dy). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&dz, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (dz). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&xmin, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (xmin). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&ymin, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (ymin). Exiting.\n");
      exit (4);
    }
  ret_it = fread (&zmin, sizeof (double), 1, fin);
  if (ret_it == 0)
    {
      printf ("Failure in reading data (zmin). Exiting.\n");
      exit (4);
    }

  data = (double *)malloc (sizeof (double) * (19 + 3 * np + 10));

  if (data == NULL)
    {
      printf ("Sorry, but I cannot allocate the data array to temporary store the input data...\n");
      exit (3);
    }

  fprintf (fout, "number of events: %12ld  \n", nevents);
  fprintf (fout, "time: %12.3e  \n", time);
  fprintf (fout, "number of particle species: %2d  \n", np);
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
	      bp=0;
              ret_it = fread (data, sizeof (double), 19 + 3 * np + 10, fin);
              if (ret_it == 0)
                {
                  printf ("Failure in reading data. Exiting.\n");
                  exit (4);
                }
              fprintf (fout, "position: %7.3e  %7.3e  %7.3e  \n", xmin + (i + 0.5) * dx, ymin + (j + 0.5) * dy,
                       zmin + (k + 0.5) * dz);
              fprintf (fout, "four velocity components u0-u4: %14.9e  %14.9e  %14.9e  %14.9e\n",
		       data[bp], data[bp+1], data[bp+2], data[bp+3]);
	      fprintf (fout, "three velocity components vx,vy vz: %14.9e  %14.9e  %14.9e\n", data[bp+1]/data[bp],
			     data[bp+2]/data[bp], data[bp+3]/data[bp]);
	      bp+=4;
	      fprintf (fout, "Ib_diffusion components: %14.9e  %14.9e  %14.9e  %14.9e\n",
		       data[bp], data[bp+1], data[bp+2], data[bp+3]);
	      bp+=4;
	      fprintf (fout, "Ic_diffusion components: %14.9e  %14.9e  %14.9e  %14.9e\n",
		       data[bp], data[bp+1], data[bp+2], data[bp+3]);
	      bp+=4;
	      fprintf (fout, "Is_diffusion components: %14.9e  %14.9e  %14.9e  %14.9e\n",
		       data[bp], data[bp+1], data[bp+2], data[bp+3]);
	      bp+=4;
	      fprintf (fout, "Particle ID, number of particles, density, energy density:\n");
              for (p = 0; p < np; p++)
                {
		   fprintf (fout, "%d  %14.9e  %14.9e  %14.9e\n", p, data[bp+p], data[bp+p+1], data[bp+p+2]);
                   bp+=3;		   
                }
	      fprintf (fout, "T00,T01,T02,T03,T11,T12,T13,T22,T23,T33:\n");
              for (p = 0; p < 10; p++)
                {
		   fprintf (fout, "%14.9e  ", data[bp+p]);
                }
	      fprintf (fout, "\n");
            }
        }
    }

  free (data);
  fclose (fout);
  fclose (fin);
  return 0;
}
