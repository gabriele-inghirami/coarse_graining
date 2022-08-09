#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>


void help()
{
printf("Syntax: ./to_landau_test.exe\n");
}

int main(int argc, char* argv[])
{

gsl_matrix *Tmunu, *lambda_mat, *Tmunu_Landau, *tmp_mat;
gsl_matrix_complex *result;
gsl_eigen_nonsymmv_workspace *ws;
gsl_vector_complex *eval;
gsl_complex eigenval, num;
gsl_vector_complex_view eigenvec;
double v[4], T_se[10];
double norm, norm2, ev;
double glf1;

Tmunu=gsl_matrix_alloc(4, 4);
if(Tmunu==NULL) {
  printf("Unable to allocate the gsl_matrix Tmunu 2D (4,4) array. I am forced to quit.\n");
  exit(4);
} 

result=gsl_matrix_complex_alloc(4, 4);
if(result==NULL) {
  printf("Unable to allocate the gsl_matrix_complex result 2D (4,4) array. I am forced to quit.\n");
  exit(4);
} 

ws=gsl_eigen_nonsymmv_alloc(4);
if(ws==NULL) {
  printf("Unable to allocate the gsl_eigen_nonsymmv ws array. I am forced to quit.\n");
  exit(4);
} 

eval=gsl_vector_complex_alloc(4);
if(eval==NULL) {
  printf("Unable to allocate the gsl_vector_complex eval array. I am forced to quit.\n");
  exit(4);
} 

if(argc==4) {
     lambda_mat=gsl_matrix_alloc(4, 4);
     if(lambda_mat==NULL) {
          printf("Unable to allocate the gsl_matrix lambda_mat 2D (4,4) array. I am forced to quit.\n");
          exit(4);
     }
     Tmunu_Landau=gsl_matrix_alloc(4, 4);
     if(Tmunu_Landau==NULL) {
          printf("Unable to allocate the gsl_matrix lambda_mat 2D (4,4) array. I am forced to quit.\n");
          exit(4);
     }
     tmp_mat=gsl_matrix_alloc(4, 4);
     if(tmp_mat==NULL) {
          printf("Unable to allocate the gsl_matrix tmp_mat 2D (4,4) array. I am forced to quit.\n");
          exit(4);
     }
} 


fread(&nevents,sizeof(long int),1,fin);
fread(&time,sizeof(double),1,fin);
fread(&np,sizeof(int),1,fin);
fread(&nx,sizeof(int),1,fin);
fread(&ny,sizeof(int),1,fin);
fread(&nz,sizeof(int),1,fin);
fread(&dx,sizeof(double),1,fin);
fread(&dy,sizeof(double),1,fin);
fread(&dz,sizeof(double),1,fin);
fread(&xmin,sizeof(double),1,fin);
fread(&ymin,sizeof(double),1,fin);
fread(&zmin,sizeof(double),1,fin);

Tp=(double *)malloc(nx*ny*nz*np*10*sizeof(double));
if(Tp==NULL)
{
 printf("Sorry, but it is not possible to allocate the Tmunu array. I am forced to quit.\n");
 exit(4);
}
Jp=(double *)malloc(nx*ny*nz*np*4*sizeof(double));
if(Jp==NULL)
{
 printf("Sorry, but it is not possible to allocate the Jp array. I am forced to quit.\n");
 exit(4);
}
Jb=(double *)malloc(nx*ny*nz*4*sizeof(double));
if(Jb==NULL)
{
 printf("Sorry, but it is not possible to allocate the Jb array. I am forced to quit.\n");
 exit(4);
}
Jc=(double *)malloc(nx*ny*nz*4*sizeof(double));
if(Jc==NULL)
{
 printf("Sorry, but it is not possible to allocate the Jc array. I am forced to quit.\n");
 exit(4);
}
Js=(double *)malloc(nx*ny*nz*4*sizeof(double));
if(Js==NULL)
{
 printf("Sorry, but it is not possible to allocate the Js array. I am forced to quit.\n");
 exit(4);
}
Jt=(double *)malloc(nx*ny*nz*4*sizeof(double));
if(Jt==NULL)
{
 printf("Sorry, but it is not possible to allocate the Jt array. I am forced to quit.\n");
 exit(4);
}
Pnum=(long int *)malloc(nx*ny*nz*np*sizeof(long int));
if(Pnum==NULL)
{
 printf("Sorry, but it is not possible to allocate the Pnum array inside main. I am forced to quit.\n");
 exit(4);
}
empty_arr=(double *)calloc((15+3*np),sizeof(double));
if(empty_arr==NULL)
{
 printf("Sorry, but it is not possible to allocate the empty_arr array inside main. I am forced to quit.\n");
 exit(4);
}

  
gsl_matrix_set(Tmunu, 0, 0, T_se[T00]);
gsl_matrix_set(Tmunu, 0, 1, -T_se[T01]);
gsl_matrix_set(Tmunu, 0, 2, -T_se[T02]);
gsl_matrix_set(Tmunu, 0, 3, -T_se[T03]);
gsl_matrix_set(Tmunu, 1, 0, T_se[T01]);
gsl_matrix_set(Tmunu, 1, 1, -T_se[T11]);
gsl_matrix_set(Tmunu, 1, 2, -T_se[T12]);
gsl_matrix_set(Tmunu, 1, 3, -T_se[T13]);
gsl_matrix_set(Tmunu, 2, 0, T_se[T02]);
gsl_matrix_set(Tmunu, 2, 1, -T_se[T12]);
gsl_matrix_set(Tmunu, 2, 2, -T_se[T22]);
gsl_matrix_set(Tmunu, 2, 3, -T_se[T23]);
gsl_matrix_set(Tmunu, 3, 0, T_se[T03]);
gsl_matrix_set(Tmunu, 3, 1, -T_se[T13]);
gsl_matrix_set(Tmunu, 3, 2, -T_se[T23]);
gsl_matrix_set(Tmunu, 3, 3, -T_se[T33]);

gsl_eigen_nonsymmv (Tmunu, eval, result,ws);

stop=0;
success=0;
for(mm=0;mm<4;mm++)
{
    go_on=0;
    eigenval=gsl_vector_complex_get(eval,mm);
    ev=GSL_REAL(eigenval);
    if((ev<=0) || (GSL_IMAG(eigenval)!=0)) continue;
    eigenvec=gsl_matrix_complex_column (result,mm);
    for(nn=0;nn<4;nn++)
	{
         num=gsl_vector_complex_get(&eigenvec.vector,nn);
         if(GSL_IMAG(num)!=0)
         {
             go_on=1;
             break;
         }
         v[nn]=GSL_REAL(num);
    }
    if((v[0]<=0) || (go_on == 1)) continue;
    norm2=v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3];
    if(norm2 <= 0) continue;
    norm=sqrt(norm2);
    for(nn=0;nn<4;nn++) u4[nn]=v[nn]/norm;
    success=1;
    break;
}
if (success!=1) {
    printf("Failure...\n");
    return 1;
}
else {
    printf("Success:\n");
    for(nn=0;nn<4;nn++) printf("%lf\n",u4[nn]);
    return 0;
}
}
