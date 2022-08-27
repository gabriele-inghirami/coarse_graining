#define TLOC 10 * (p + np * (k + nz * (j + ny * (i + h * nx))))
#define JPL 4 * (p + np * (k + nz * (j + ny * (i + h * nx))))
#define JBL 4 * (k + nz * (j + ny * (i + h * nx)))
#define PNLOC (p + np * (k + nz * (j + ny * (i + (long)h * nx))))
#define TLOC2 10 * (p + np * (k2 + nz2 * (j2 + ny2 * (i2 + h * nx2))))
#define JPL2 4 * (p + np * (k2 + nz2 * (j2 + ny2 * (i2 + h * nx2))))
#define JBL2 4 * (k2 + nz2 * (j2 + ny2 * (i2 + h * nx2)))
#define PNLOC2 (p + np * (k2 + nz2 * (j2 + ny2 * (i2 + (long)h * nx2))))

const int nx2 = 7;
const int ny2 = 7;
const int nz2 = 7;
const int define_grid_corner; // if different from 0 the program uses xmin2, ymin2, zmin2 as grid corner
const double new_xmin = -1;
const double new_ymin = -1;
const double new_zmin = -1;

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
}; /* The flattened indexes of the energy momentum tensor components */
enum
{
  J0 = 0,
  J1,
  J2,
  J3
}; /* The flattened indexes of four momentum components */
const int T10 = T01;
const int T20 = T02;
const int T30 = T03;
const int T21 = T12;
const int T31 = T13;
const int T32 = T23;

#include <stdio.h>
#include <stdlib.h>

void
help ()
{
  printf ("Syntax: ./resize.exe <Tmunu_inputfile> <Tmunu_outputfile>\n\n");
}

int
main (int argc, char *argv[])
{

  FILE *fin, *fTp;
  double *Tp, *Jb, *Jp, *Jc, *Js, *Jt, *Tp2, *Jb2, *Jp2, *Jc2, *Js2, *Jt2;
  long int *Pnum, *Pnum2;
  long int nevents;
  double time;
  int i, j, k, p, l, i2, j2, k2;
  int nx, ny, nz, np;
  double dx, dy, dz, xmin, ymin, zmin, xmin2, ymin2, zmin2;
  const int h = 0;

  int i_start, i_end, j_start, j_end, k_start, k_end;

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

  fTp = fopen (argv[2], "wb");
  if (fTp == NULL)
    {
      printf ("Sorry, I was unable to create the output file %s\n", argv[2]);
      exit (2);
    }

  fread (&nevents, sizeof (long int), 1, fin);
  fread (&time, sizeof (double), 1, fin);
  fread (&np, sizeof (int), 1, fin);
  fread (&nx, sizeof (int), 1, fin);
  fread (&ny, sizeof (int), 1, fin);
  fread (&nz, sizeof (int), 1, fin);
  fread (&dx, sizeof (double), 1, fin);
  fread (&dy, sizeof (double), 1, fin);
  fread (&dz, sizeof (double), 1, fin);
  fread (&xmin, sizeof (double), 1, fin);
  fread (&ymin, sizeof (double), 1, fin);
  fread (&zmin, sizeof (double), 1, fin);

  Tp = (double *)malloc (nx * ny * nz * np * 10 * sizeof (double));
  if (Tp == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Tmunu array. I am forced to quit.\n");
      exit (4);
    }
  Jp = (double *)malloc (nx * ny * nz * np * 4 * sizeof (double));
  if (Jp == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jp array. I am forced to quit.\n");
      exit (4);
    }
  Jb = (double *)malloc (nx * ny * nz * 4 * sizeof (double));
  if (Jb == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jb array. I am forced to quit.\n");
      exit (4);
    }
  Jc = (double *)malloc (nx * ny * nz * 4 * sizeof (double));
  if (Jc == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jc array. I am forced to quit.\n");
      exit (4);
    }
  Js = (double *)malloc (nx * ny * nz * 4 * sizeof (double));
  if (Js == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Js array. I am forced to quit.\n");
      exit (4);
    }
  Jt = (double *)malloc (nx * ny * nz * 4 * sizeof (double));
  if (Jt == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jt array. I am forced to quit.\n");
      exit (4);
    }
  Pnum = (long int *)malloc (nx * ny * nz * np * sizeof (long int));
  if (Pnum == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Pnum array inside main. I am forced to quit.\n");
      exit (4);
    }

  Tp2 = (double *)malloc (nx2 * ny2 * nz2 * np * 10 * sizeof (double));
  if (Tp2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Tmunu2 array. I am forced to quit.\n");
      exit (4);
    }
  Jp2 = (double *)malloc (nx2 * ny2 * nz2 * np * 4 * sizeof (double));
  if (Jp2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jp2 array. I am forced to quit.\n");
      exit (4);
    }
  Jb2 = (double *)malloc (nx2 * ny2 * nz2 * 4 * sizeof (double));
  if (Jb2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jb2 array. I am forced to quit.\n");
      exit (4);
    }
  Jc2 = (double *)malloc (nx2 * ny2 * nz2 * 4 * sizeof (double));
  if (Jc2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jc2 array. I am forced to quit.\n");
      exit (4);
    }
  Js2 = (double *)malloc (nx2 * ny2 * nz2 * 4 * sizeof (double));
  if (Js2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Js2 array. I am forced to quit.\n");
      exit (4);
    }
  Jt2 = (double *)malloc (nx2 * ny2 * nz2 * 4 * sizeof (double));
  if (Jt2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Jt2 array. I am forced to quit.\n");
      exit (4);
    }
  Pnum2 = (long int *)malloc (nx2 * ny2 * nz2 * np * sizeof (long int));
  if (Pnum2 == NULL)
    {
      printf ("Sorry, but it is not possible to allocate the Pnum array inside main. I am forced to quit.\n");
      exit (4);
    }

  fread (Tp, sizeof (double), nx * ny * nz * np * 10, fin);
  fread (Jp, sizeof (double), nx * ny * nz * np * 4, fin);
  fread (Jb, sizeof (double), nx * ny * nz * 4, fin);
  fread (Jc, sizeof (double), nx * ny * nz * 4, fin);
  fread (Js, sizeof (double), nx * ny * nz * 4, fin);
  fread (Jt, sizeof (double), nx * ny * nz * 4, fin); // currently Jt is read, but not used
  fread (Pnum, sizeof (long int), nx * ny * nz * np, fin);
  printf ("%s read.\n", argv[1]);
  fclose (fin);

  // now we start to write in the output file

  fwrite (&nevents, sizeof (long int), 1, fTp);
  fwrite (&time, sizeof (double), 1, fTp);
  fwrite (&np, sizeof (int), 1, fTp);
  fwrite (&nx2, sizeof (int), 1, fTp);
  fwrite (&ny2, sizeof (int), 1, fTp);
  fwrite (&nz2, sizeof (int), 1, fTp);
  fwrite (&dx, sizeof (double), 1, fTp);
  fwrite (&dy, sizeof (double), 1, fTp);
  fwrite (&dz, sizeof (double), 1, fTp);

  if (define_grid_corner != 0)
    {
      xmin2 = new_xmin;
      ymin2 = new_ymin;
      zmin2 = new_zmin;
    }
  else
    {
      xmin2 = -nx * dx / 2;
      ymin2 = -ny * dy / 2;
      zmin2 = -nz * dz / 2;
    }

  for (i = 0; i < nx; i++)
    {
      if (xmin + (i + 0.5) * dx > xmin2)
        {
          i_start = i;
          i_end = i_start + nx2;
          xmin2 = xmin + i * dx;
          break;
        }
    }

  for (i = 0; i < ny; i++)
    {
      if (ymin + (i + 0.5) * dy > ymin2)
        {
          j_start = i;
          j_end = j_start + ny2;
          ymin2 = ymin + i * dy;
          break;
        }
    }

  for (i = 0; i < nz; i++)
    {
      if (zmin + (i + 0.5) * dz > zmin2)
        {
          k_start = i;
          k_end = k_start + nz2;
          zmin2 = zmin + i * dz;
          break;
        }
    }

  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (p = 0; p < np; p++)
                {
                  for (l = 0; l < 10; l++)
                    Tp2[l + TLOC2] = Tp[l + TLOC];
                }
            }
        }
    }

  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (p = 0; p < np; p++)
                {
                  for (l = 0; l < 4; l++)
                    Jp2[l + JPL2] = Jp[l + JPL];
                }
            }
        }
    }

  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (l = 0; l < 4; l++)
                Jb2[l + JBL2] = Jb[l + JBL];
            }
        }
    }
  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (l = 0; l < 4; l++)
                Jc2[l + JBL2] = Jc[l + JBL];
            }
        }
    }
  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (l = 0; l < 4; l++)
                Js2[l + JBL2] = Js[l + JBL];
            }
        }
    }
  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (l = 0; l < 4; l++)
                Jt2[l + JBL2] = Jt[l + JBL];
            }
        }
    }

  for (i = i_start; i < i_end; i++)
    {
      i2 = i - i_start; // offset along i
      for (j = j_start; j < j_end; j++)
        {
          j2 = j - j_start; // offset along j
          for (k = k_start; k < k_end; k++)
            {
              k2 = k - k_start; // offset along k
              for (p = 0; p < np; p++)
                Pnum2[PNLOC2] = Pnum[PNLOC];
            }
        }
    }
  fwrite (&xmin2, sizeof (double), 1, fTp);
  fwrite (&ymin2, sizeof (double), 1, fTp);
  fwrite (&zmin2, sizeof (double), 1, fTp);

  fwrite (Tp2, sizeof (double), nx2 * ny2 * nz2 * np * 10, fTp);
  fwrite (Jp2, sizeof (double), nx2 * ny2 * nz2 * np * 4, fTp);
  fwrite (Jb2, sizeof (double), nx2 * ny2 * nz2 * 4, fTp);
  fwrite (Jc2, sizeof (double), nx2 * ny2 * nz2 * 4, fTp);
  fwrite (Js2, sizeof (double), nx2 * ny2 * nz2 * 4, fTp);
  fwrite (Jt2, sizeof (double), nx2 * ny2 * nz2 * 4, fTp);
  fwrite (Pnum2, sizeof (long int), nx2 * ny2 * nz2 * np, fTp);

  printf ("Data of the resized grid saved in file %s.\n", argv[2]);
  fclose (fTp);

  free (Tp);
  free (Jb);
  free (Jp);
  free (Jc);
  free (Js);
  free (Jt);
  free (Pnum);
  free (Tp2);
  free (Jb2);
  free (Jp2);
  free (Jc2);
  free (Js2);
  free (Jt2);
  free (Pnum2);
  return 0;
}
