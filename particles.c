#include "definitions.h"

/** @file particles.c
 *
 *   @brief this file contains the functions dealing with the particles to include in the computation:
 * get_particle_index and associate_particle_array_index
 */

#ifdef URQMD
// given the UrQMD particle itpye and iso3, it return the particle index
int
get_particle_index (int itype, int iso3)
{
  int result;
  // Warning: this index has to be manually updated and, every time it is modified, the function
  // associate_particle_array_index must be modified, too.
  switch (itype)
    {
    case (101):
      if (iso3 == -2)
        result = 0; // pion -
      else if (iso3 == 2)
        result = 1; // pion +
      else
        result = 2; // pion 0
      break;        // the break statemets here are not really useful, but not even dangerous.
    case (106):
      if (iso3 == -1)
        result = 3; // Kaon 0
      else
        result = 4; // Kaon +
      break;
    case (-106):
      if (iso3 == -1)
        result = 5; // Kaon -
      else
        result = 6; // Kaon 0 bar
      break;
    case (1):
      if (iso3 == -1)
        result = 7; // Neutron
      else
        result = 8; // Proton
      break;
    case (-1):
      if (iso3 == -1)
        result = 9; // Anti-proton
      else
        result = 10; // Anti-neutron
      break;
    case (102):
      result = 11; // eta
      break;
    case (103):
      result = 12; // omega meson
      break;
    case (107):
      result = 13; // eta'
      break;
    case (109):
      result = 14; // phi
      break;
    case (27):
      result = 15; // Lambda1116
      break;
    case (-27):
      result = 16; // Anti-Lambda1116
      break;
    case (40):
      // Sigma1192
      if (iso3 == -2)
        result = 17; // Sigma -
      else if (iso3 == 2)
        result = 18; // Sigma +
      else
        result = 19; // Sigma 0
      break;
    case (-40):
      // Anti-Sigma1192
      if (iso3 == -2)
        result = 20; // Anti-Sigma (charge -)
      else if (iso3 == 2)
        result = 21; // Anti-Sigma (charge +)
      else
        result = 22; // Anti-Sigma 0
      break;
    case (49):
      // Xi1317
      if (iso3 == -1)
        result = 23; // Xi -
      else
        result = 24; // Xi 0
      break;
    case (-49):
      // Anti-Xi1317
      if (iso3 == -1)
        result = 25; // Anti-Xi 0
      else
        result = 26; // Xi +
      break;
    case (29):
      result = 27; // Lambda1520
      break;
    case (-29):
      result = 28; // Anti-Lambda1520
      break;
    case (50):
      // Xi1530
      if (iso3 == -1)
        result = 29; // Xi -
      else
        result = 30; // Xi 0
      break;         // Xi1530
    case (-50):
      // Anti-Xi1530
      if (iso3 == -1)
        result = 31; // Anti-Xi 0
      else
        result = 32; // Xi +
      break;
    case (55):
      result = 33; // Omega1672
      break;
    case (-55):
      result = 34; // Anti-Omega1672
      break;
    // index for other particles
    default:
      result = NP; // all other particles
    }
  if (result > NP)
    result = NP;
  return result;
}

#elif defined(SMASH)

/**
 * plist_unsorted (SMASH)
 * The array containing the characteristics of the particles managed by SMASH.
 * It has a fixed dimension that should be sufficient even for the upcoming versions of SMASH.
 * This array is later ordered by growing pdg_id and the results are store in the array plist
 */
pinfo plist_unsorted[MAX_SMASH_HADRON_SPECIES];

/**
 * plist (SMASH)
 * This is a pointer the first element of an array containing the characteristics of the particles managed by SMASH
 * This array is ordered by growing pdg_id and dinamically allocated at runtime
 */
pinfo *plist;

int n_hadrons_smash;

// given the SMASH particle pdg_id, it return the particle index
int
get_particle_index (int pdg_id)
{
  int result;
  // Warning: this index has to be manually updated and, every time it is modified, the function
  // associate_particle_array_index must be modified, too.
  switch (pdg_id)
    {
    case (-211):
      result = 0; // pion -
      break;
    case (211):
      result = 1; // pion +
      break;
    case (111):
      result = 2; // pion 0
      break;
    case (311):
      result = 3; // Kaon 0
      break;
    case (321):
      result = 4; // Kaon +
      break;
    case (-321):
      result = 5; // Kaon -
      break;
    case (-311):
      result = 6; // Kaon 0 bar
      break;
    case (2112):
      result = 7; // Neutron
      break;
    case (2212):
      result = 8; // Proton
      break;
    case (-2212):
      result = 9; // Anti-proton
      break;
    case (-2112):
      result = 10; // Anti-neutron
      break;
    case (221):
      result = 11; // eta
      break;
    case (223):
      result = 12; // omega meson
      break;
    case (331):
      result = 13; // eta'
      break;
    case (333):
      result = 14; // phi
      break;
    case (3122):
      result = 15; // Lambda1116
      break;
    case (-3122):
      result = 16; // Anti-Lambda1116
      break;
    case (3112):
      // Sigma1192
      result = 17; // Sigma -
      break;
    case (3222):
      result = 18; // Sigma +
      break;
    case (3212):
      result = 19; // Sigma 0
      break;
    case (-3222):
      // Anti-Sigma1192
      result = 20; // Anti-Sigma (charge -)
      break;
    case (-3112):
      result = 21; // Anti-Sigma (charge +)
      break;
    case (-3212):
      result = 22; // Anti-Sigma 0
      break;
    case (3312):
      // Xi1317
      result = 23; // Xi -
      break;
    case (3322):
      result = 24; // Xi 0
      break;
    case (-3322):
      // Anti-Xi1317
      result = 25; // Anti-Xi 0
      break;
    case (-3312):
      result = 26; // Xi +
      break;
    case (3124):
      result = 27; // Lambda1520
      break;
    case (-3124):
      result = 28; // Anti-Lambda1520
      break;
    case (3314):
      // Xi1530
      result = 29; // Xi -
      break;
    case (3324):
      result = 30; // Xi 0
      break;       // Xi1530
    case (-3324):
      // Anti-Xi1530
      result = 31; // Anti-Xi 0
      break;
    case (-3314):
      result = 32; // Xi +
      break;
    case (3334):
      result = 33; // Omega1672
      break;
    case (-3334):
      result = 34; // Anti-Omega1672
      break;
    // index for other particles
    default:
      result = NP; // all other particles
    }
  if (result > NP)
    result = NP;
  return result;
}

#endif

// it returns the name of the particle given its array index
const char *
associate_particle_array_index (int index)
{
  switch (index)
    {
    case (0):
      return "Pion- * UrQMD itype 101, 2iso3 -2, charge -1 * PDG ID -211";
      break;
    case (1):
      return "Pion+ * UrQMD itype 101, 2iso3 2, charge +1 * PDG ID 211";
      break;
    case (2):
      return "Pion0 * UrQMD itype 101, 2iso3 0, charge 0 * PDG ID 111";
      break;
    case (3):
      return "Kaon0 * UrQMD itype 106, 2iso3 -1, charge 0 * PDG ID 311";
      break;
    case (4):
      return "Kaon+ * UrQMD itype 106, 2iso3 1, charge +1 * PDG ID 321";
      break;
    case (5):
      return "Kaon- * UrQMD itype -106, 2iso3 -1, charge -1 * PDG ID -321";
      break;
    case (6):
      return "Kaon0 bar * UrQMD itype -106, 2iso3 1, charge 0 * PDG ID -311";
      break;
    case (7):
      return "Neutron * UrQMD itype 1, 2iso3 -1, charge 0 * PDG ID 2112";
      break;
    case (8):
      return "Proton * UrQMD itype 1, 2iso3 1, charge +1 * PDG ID 2212";
      break;
    case (9):
      return "Anti-Proton * UrQMD itype -1, 2iso3 -1, charge -1 * PDG ID -2212";
      break;
    case (10):
      return "Anti-Neutron * UrQMD itype -1, 2iso3 1, charge 0 * PDG ID -2112";
      break;
    case (11):
      return "Eta meson * UrQMD itype 102, 2iso3 0, charge 0 * PDG ID 221";
      break;
    case (12):
      return "omega meson * UrQMD itype 103, 2iso3 0, charge 0 * PDG ID 223";
      break;
    case (13):
      return "eta' meson * UrQMD itype 107, 2iso3 0, charge 0 * PDG ID 331";
      break;
    case (14):
      return "phi meson * UrQMD itype 109, 2iso3 0, charge 0 * PDG ID 333";
      break;
    case (15):
      return "Lambda1116 * UrQMD itype 27, 2iso3 0, charge 0 * PDG ID 3122";
      break;
    case (16):
      return "Anti-Lambda1116 * UrQMD itype -27, 2iso3 0, charge 0 * PDG ID -3122";
      break;
    case (17):
      return "Sigma1192- * UrQMD itype 40, 2iso3 -2, charge -1, dds * PDG ID 3112";
      break;
    case (18):
      return "Sigma1192+ * UrQMD itype 40, 2iso3 2, charge 1, uus * PDG ID 3222";
      break;
    case (19):
      return "Sigma1192 * UrQMD itype 40, 2iso3 0, charge 0, uds * PDG ID 3212";
      break;
    case (20):
      return "Anti-Sigma1192- * UrQMD itype -40, 2iso3 -2, charge -1, u_bar u_bar s_bar * PDG ID -3222";
      break;
    case (21):
      return "Anti-Sigma1192+ * UrQMD itype -40, 2iso3 2, charge +1, d_bar d_bar s_bar * PDG ID -3112";
      break;
    case (22):
      return "Anti-Sigma1192 * UrQMD itype -40, 2iso3 0, charge 0, u_bar d_bar s_bar * PDG ID -3212";
      break;
    case (23):
      return "Xi1317- * UrQMD itype 49, 2iso3 -1, charge -1 * PDG ID 3312";
      break;
    case (24):
      return "Xi1317 0 * UrQMD itype 49, 2iso3 1, charge 0 * PDG ID 3322";
      break;
    case (25):
      return "Anti-Xi1317 0 * UrQMD itype -49, 2iso3 -1, charge 0 * PDG ID -3322";
      break;
    case (26):
      return "Xi1317+ * UrQMD itype -49, 2iso3 1, charge +1 * PDG ID -3312";
      break;
    case (27):
      return "Lambda1520 * UrQMD itype 29, 2iso3 0, charge 0 * PDG ID 3124";
      break;
    case (28):
      return "Anti-Lambda1520 * UrQMD itype -29, 2iso3 0, charge 0 * PDG ID -3124";
      break;
    case (29):
      return "Xi1530- * UrQMD itype 50, 2iso3 -1, charge -1 * PDG ID 3314";
      break;
    case (30):
      return "Xi1530 0 * UrQMD itype 50, 2iso3 1, charge 0 * PDG ID 3324";
      break;
    case (31):
      return "Anti-Xi1530 0 * UrQMD itype -50, 2iso3 -1, charge 0 * PDG ID -3324";
      break;
    case (32):
      return "Xi1530+ * UrQMD itype -50, 2iso3 1, charge +1 * PDG ID -3314";
      break;
    case (33):
      return "Omega1672 * UrQMD itype 55, 2iso3 0, charge -1 * PDG ID 3334";
      break;
    case (34):
      return "Anti-Omega1672 * UrQMD itype -55, 2iso3 0, charge +1 * PDG ID -3334";
      break;
    default:
      return "All other particles";
    }
}

#ifdef URQMD
// it returns the strangeness given the itype
int
get_strangeness (int itype)
{

  int at;

  at = abs (itype);

  if (at < 27)
    {
      return 0;
    }
  else if (at < 49)
    {
      if (itype > 0)
        return -1;
      return 1;
    }
  else if (at < 55)
    {
      if (itype > 0)
        return -2;
      return 2;
    }
  else if (at == 55)
    {
      if (itype > 0)
        return -3;
      return 3;
    }
  else if (at < 106)
    {
      return 0;
    }
  else if (at == 106)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at == 107)
    {
      return 0;
    }
  else if (at == 108)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at == 109)
    {
      return 0;
    }
  else if (at == 110)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at < 113)
    {
      return 0;
    }
  else if (at == 113)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at < 117)
    {
      return 0;
    }
  else if (at == 117)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at < 121)
    {
      return 0;
    }
  else if (at == 121)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at < 125)
    {
      return 0;
    }
  else if (at == 125)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at < 129)
    {
      return 0;
    }
  else if (at == 129)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else if (at < 138)
    {
      return 0;
    }
  else if (at < 140)
    {
      if (itype > 0)
        return 1;
      return -1;
    }
  else
    {
      printf ("Unidentified particle with itype %d. Assigning 0 strangeness.\n", itype);
      return 0;
    }
}

#elif defined(SMASH)
void
get_had_prop (char *id_string, int id_string_len, int *s, int *B, int *Jtot)
{
  int q1, q2, q3, q3_tmp, end;
  int strangeness[] = { 0, 0, -1, 0, 0 };
  double c_tmp;
  end = id_string_len - 1;
  if (end == 1) // they are leptons
    {
      *s = 0;
      *B = 0;
      *Jtot = 2;
      return;
    }

  // printf("Lenght: %d\n",end);
  // printf("Id: %s\n",id_string);
  // printf("%c\n",id_string[end]);
  *Jtot = (int)id_string[end] - (int)'0';
  q1 = (int)id_string[end - 1] - (int)'0';
  q2 = (int)id_string[end - 2] - (int)'0';
  if (end > 2)
    {
      q3_tmp = (int)id_string[end - 3] - (int)'0';
      if (q3_tmp > 0)
        {
          q3 = q3_tmp;
          *B = 1; // 3 quarks: it is a baryon
          if (*Jtot == 9)
            *Jtot = 10; // this is an expection for N(2200), N(2250) and Lambda(2450) with J=9/2, so that 2J+1=10, but
                        // in their ID the last number is 9
          if ((*Jtot % 2) != 0)
            { // it must be a fermion, 2J+1 must be even and, if divided by 2 there is a reminder, it is odd...
              printf ("Sorry, but there is something wrong. In particles.txt I have found a baryon with integer "
                      "spin... I quit.\n");
              exit (6);
            }
        }
      else
        {
          q3 = 0;
          *B = 0; // 2 quarks: it is a meson
          if ((*Jtot % 2) == 0)
            { // it must be a boson, 2J+1 must be odd
              printf ("Sorry, but there is something wrong. In particles.txt I found a meson with half-integer "
                      "spin... I quit.\n");
              exit (6);
            }
        }
    }
  else
    q3 = 0;

  // we shift the indexes by 1
  *s = strangeness[q1 - 1] + strangeness[q2 - 1];
  if (q3 > 0)
    *s = *s + strangeness[q3 - 1];
  else
    *s = -*s; // if q3==0 is a meson and in that case the s are actually anti-s
  if ((q3 == 0) && (q1 == q2))
    *s = 0; // it is a s-antis state
}

int
fill (pinfo *plist)
{
  int n = 0;
  int i;
  wchar_t buf[200];
  wchar_t tmp[16];
  FILE *finp;
  float float_bb;
  char par;
  char pdg_id_stringA[9], pdg_id_stringB[9], pdg_id_stringC[9], pdg_id_stringD[9];
  const char zero_string[9] = { '0', '0', '0', '0', '0', '0', '0', '0', '\n' };
  int tmp_pdg_idB, tmp_pdg_idC, tmp_pdg_idD;

  setlocale (LC_ALL, "en_US.utf8");
  finp = fopen ("particles.txt", "r");
  if (finp == NULL)
    {
      printf ("Sorry, but it is not possible to open particles.txt, containing the data of the hadron species managed "
              "by SMASH... I quit.\n\n");
      exit (5);
    }

  while (fgetws (buf, 200, finp) != NULL)
    {
      if ((buf[0] == '#') || (buf[0] == ' ') || (buf[0] == '\n'))
        {
          continue;
        }
      else
        {
          swscanf (buf, L"%ls %lf %e %c %s %s %s %s", &plist[n].name, &plist[n].mass, &float_bb, &par, &pdg_id_stringA,
                   &pdg_id_stringB, &pdg_id_stringC, &pdg_id_stringD);
          plist[n].pdg_id = atoi (pdg_id_stringA);
          get_had_prop (pdg_id_stringA, strlen (pdg_id_stringA), &plist[n].strangeness, &plist[n].baryon_num,
                        &plist[n].spin);
          // printf("%-3d   %-11ls  %12.4lf %10d   B: %1d   s: %2d   2J+1: %2d\n",n, plist[n].name, plist[n].mass,
          // plist[n].pdg_id, plist[n].baryon_num, plist[n].strangeness, plist[n].spin);
          n = n + 1;

          if ((pdg_id_stringB[0] == '#') || (pdg_id_stringB[0] == ' ') || (pdg_id_stringB[0] == '\0')
              || (pdg_id_stringB[0] == '\n'))
            continue;
          tmp_pdg_idB = atoi (pdg_id_stringB);
          if (tmp_pdg_idB != 0)
            {
              wcsncpy ((int *)&plist[n].name, (int *)&plist[n - 1].name, wcslen (plist[n - 1].name));
              plist[n].mass = plist[n - 1].mass;
              plist[n].pdg_id = tmp_pdg_idB;
              get_had_prop (pdg_id_stringB, strlen (pdg_id_stringB), &plist[n].strangeness, &plist[n].baryon_num,
                            &plist[n].spin);
              // printf("%-3d   %-11ls  %12.4lf %10d   B: %1d   s: %2d   2J+1: %2d\n",n, plist[n].name, plist[n].mass,
              // plist[n].pdg_id, plist[n].baryon_num, plist[n].strangeness, plist[n].spin);
              n = n + 1;
              strncpy (pdg_id_stringB, zero_string, sizeof (pdg_id_stringB));
            }

          if ((pdg_id_stringC[0] == '#') || (pdg_id_stringC[0] == ' ') || (pdg_id_stringC[0] == '\0')
              || (pdg_id_stringC[0] == '\n'))
            continue;
          tmp_pdg_idC = atoi (pdg_id_stringC);
          if (tmp_pdg_idC != 0)
            {
              wcsncpy ((int *)&plist[n].name, (int *)&plist[n - 1].name, wcslen (plist[n - 1].name));
              plist[n].mass = plist[n - 1].mass;
              plist[n].pdg_id = tmp_pdg_idC;
              get_had_prop (pdg_id_stringC, strlen (pdg_id_stringC), &plist[n].strangeness, &plist[n].baryon_num,
                            &plist[n].spin);
              // printf("%-3d   %-11ls  %12.4lf %10d   B: %1d   s: %2d   2J+1: %2d\n",n, plist[n].name, plist[n].mass,
              // plist[n].pdg_id, plist[n].baryon_num, plist[n].strangeness, plist[n].spin);
              n = n + 1;
              strncpy (pdg_id_stringC, zero_string, sizeof (pdg_id_stringC));
            }

          if ((pdg_id_stringD[0] == '#') || (pdg_id_stringD[0] == ' ') || (pdg_id_stringD[0] == '\0')
              || (pdg_id_stringD[0] == '\n'))
            continue;
          tmp_pdg_idD = atoi (pdg_id_stringD);
          if (tmp_pdg_idD != 0)
            {
              wcsncpy ((int *)&plist[n].name, (int *)&plist[n - 1].name, wcslen (plist[n - 1].name));
              plist[n].mass = plist[n - 1].mass;
              plist[n].pdg_id = tmp_pdg_idD;
              get_had_prop (pdg_id_stringD, strlen (pdg_id_stringD), &plist[n].strangeness, &plist[n].baryon_num,
                            &plist[n].spin);
              // printf("%-3d   %-11ls  %12.4lf %10d   B: %1d   s: %2d   2J+1: %2d\n",n, plist[n].name, plist[n].mass,
              // plist[n].pdg_id, plist[n].baryon_num, plist[n].strangeness, plist[n].spin);
              n = n + 1;
              strncpy (pdg_id_stringD, zero_string, sizeof (pdg_id_stringD));
            }
        }
    }

  fclose (finp);
  return n;
}

int
compare_for_sorting (const void *a, const void *b)
{
  const pinfo *pa = (const pinfo *)a;
  const pinfo *pb = (const pinfo *)b;

  return pa->pdg_id - pb->pdg_id;
}

void
prepare_smash_hadron_array ()
{
  int i;
  n_hadrons_smash = fill (plist_unsorted);
  plist = (pinfo *)malloc (sizeof (pinfo) * n_hadrons_smash);
  if (plist == NULL)
    {
      printf ("Unable to allocate array plist in function prepare_smash_hadron_array. I quit.\n");
      exit (4);
    }
  for (i = 0; i < n_hadrons_smash; i++)
    {
      plist[i].pdg_id = plist_unsorted[i].pdg_id;
      wcsncpy ((int *)&plist[i].name, (int *)&plist_unsorted[i].name, wcslen (plist_unsorted[i].name));
      plist[i].mass = plist_unsorted[i].mass;
      plist[i].spin = plist_unsorted[i].spin;
      plist[i].strangeness = plist_unsorted[i].strangeness;
      plist[i].baryon_num = plist_unsorted[i].baryon_num;
    }

  qsort (plist, n_hadrons_smash, sizeof (pinfo), compare_for_sorting);
}

void
get_hadron_info (const int pdg_id, int *s, int *B)
{
  int sign;
  sign = (pdg_id > 0 ? 1 : -1);
  pinfo target, *result;
  target.pdg_id = abs (pdg_id);
  result = (pinfo *)bsearch (&target, plist, n_hadrons_smash, sizeof (pinfo), compare_for_sorting);
  *s = sign * result->strangeness;
  *B = sign * result->baryon_num;
}
#endif
