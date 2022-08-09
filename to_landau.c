#define TLOC 10*(p+np*(k+nz*(j+ny*(i+h*nx))))
#define JPL 4*(p+np*(k+nz*(j+ny*(i+h*nx))))
#define JBL 4*(k+nz*(j+ny*(i+h*nx)))
#define PNLOC (p+np*(k+nz*(j+ny*(i+(long)h*nx))))

enum {T00=0,T01,T02,T03,T11,T12,T13,T22,T23,T33}; /* The flattened indexes of the energy momentum tensor components */
enum {J0=0,J1,J2,J3};  /* The flattened indexes of four momentum components */
const int T10=T01;
const int T20=T02;
const int T30=T03;
const int T21=T12;
const int T31=T13;
const int T32=T23;

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
printf("Syntax: ./to_landau.exe <Tmumu_inputfile> <density_file> [Tmunu_file]\n\n");
}

int main(int argc, char* argv[])
{

FILE *fin, *fde, *ftm;
double *Tp, *Jb, *Jp, *Jc, *Js, *Jt;
long int *Pnum;
long int nevents;
double time, xmin, ymin, zmin;
int i,j,k,p,h,l,mm,nn;
int nx, ny, nz, np;
double dx, dy, dz;
double rho_c, rho_s, rho_b, rho_t;
double u4[4];
double Ic_diffusion[4], Is_diffusion[4], Ib_diffusion[4];
double Ic_diffcheck[4], Is_diffcheck[4], Ib_diffcheck[4];
double eps, rho;
gsl_matrix *Tmunu, *lambda_mat, *Tmunu_Landau, *tmp_mat;
gsl_matrix_complex *result;
gsl_eigen_nonsymmv_workspace *ws;
gsl_vector_complex *eval;
gsl_complex eigenval, num;
gsl_vector_complex_view eigenvec;
double v[4], T_se[10];
double norm, norm2, ev;
double glf1;
int stop,go_on,success;
double cell_volume,NF;
double tmp_value;
double *empty_arr;


if((argc<3) || (argc>4))
{
  help();
  exit(1);
}

fin=fopen(argv[1],"r");
if(fin==NULL)
{
	printf("Sorry, I was unable to open the input file %s\n",argv[1]);
	exit(2);
}

fde=fopen(argv[2],"wb");
if(fde==NULL)
{
	printf("Sorry, I was unable to create the density output file %s\n",argv[2]);
	exit(2);
}

if(argc==4)
{
    ftm=fopen(argv[3],"wb");
    if(ftm==NULL)
    {
	    printf("Sorry, I was unable to create the Tmunu output file %s\n",argv[3]);
	    exit(2);
    }
}

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

  
fread(Tp,sizeof(double),nx*ny*nz*np*10,fin);
fread(Jp,sizeof(double),nx*ny*nz*np*4,fin);
fread(Jb,sizeof(double),nx*ny*nz*4,fin);
fread(Jc,sizeof(double),nx*ny*nz*4,fin);
fread(Js,sizeof(double),nx*ny*nz*4,fin);
fread(Jt,sizeof(double),nx*ny*nz*4,fin); //currently Jt is read, but not used
fread(Pnum,sizeof(long int),nx*ny*nz*np,fin);
printf("%s read.\n",argv[1]);
fclose(fin);


//now we start to write in the output file				
      
fwrite(&nevents,sizeof(long int),1,fde);
fwrite(&time,sizeof(double),1,fde);
fwrite(&np,sizeof(int),1,fde);
fwrite(&nx,sizeof(int),1,fde);
fwrite(&ny,sizeof(int),1,fde);
fwrite(&nz,sizeof(int),1,fde);
fwrite(&dx,sizeof(double),1,fde);
fwrite(&dy,sizeof(double),1,fde);
fwrite(&dz,sizeof(double),1,fde);
fwrite(&xmin,sizeof(double),1,fde);
fwrite(&ymin,sizeof(double),1,fde);
fwrite(&zmin,sizeof(double),1,fde);

cell_volume=dx*dy*dz;

NF=cell_volume*nevents;

for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
      {
	for(k=0;k<nz;k++)					
          {
            for(l=0;l<10;l++) T_se[l]=0;
            for(p=0;p<np;p++) 
              {
                for(l=0;l<10;l++) T_se[l]+=Tp[l+TLOC];
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
            if(success==1) {

	      rho_b=(Jb[J0+JBL]*u4[0]-Jb[J1+JBL]*u4[1]-Jb[J2+JBL]*u4[2]-Jb[J3+JBL]*u4[3])/NF;
	      rho_c=(Jc[J0+JBL]*u4[0]-Jc[J1+JBL]*u4[1]-Jc[J2+JBL]*u4[2]-Jc[J3+JBL]*u4[3])/NF;
	      rho_s=(Js[J0+JBL]*u4[0]-Js[J1+JBL]*u4[1]-Js[J2+JBL]*u4[2]-Js[J3+JBL]*u4[3])/NF;
              fwrite(&rho_b,sizeof(double),1,fde);
              fwrite(&rho_c,sizeof(double),1,fde);
              fwrite(&rho_s,sizeof(double),1,fde);
              Ib_diffcheck[0]=(1-u4[0]*u4[0])*Jb[J0+JBL]+u4[0]*u4[1]*Jb[J1+JBL]+u4[0]*u4[2]*Jb[J2+JBL]+u4[0]*u4[3]*Jb[J3+JBL];
              Ib_diffcheck[1]=-u4[0]*u4[1]*Jb[J0+JBL]+(1+u4[1]*u4[1])*Jb[J1+JBL]+u4[1]*u4[2]*Jb[J2+JBL]+u4[1]*u4[3]*Jb[J3+JBL];
              Ib_diffcheck[2]=-u4[0]*u4[2]*Jb[J0+JBL]+u4[2]*u4[1]*Jb[J1+JBL]+(1+u4[2]*u4[2])*Jb[J2+JBL]+u4[2]*u4[3]*Jb[J3+JBL];
              Ib_diffcheck[3]=-u4[0]*u4[3]*Jb[J0+JBL]+u4[3]*u4[1]*Jb[J1+JBL]+u4[3]*u4[2]*Jb[J2+JBL]+(1+u4[3]*u4[3])*Jb[J3+JBL];
              Ic_diffcheck[0]=(1-u4[0]*u4[0])*Jc[J0+JBL]+u4[0]*u4[1]*Jc[J1+JBL]+u4[0]*u4[2]*Jc[J2+JBL]+u4[0]*u4[3]*Jc[J3+JBL];
              Ic_diffcheck[1]=-u4[0]*u4[1]*Jc[J0+JBL]+(1+u4[1]*u4[1])*Jc[J1+JBL]+u4[1]*u4[2]*Jc[J2+JBL]+u4[1]*u4[3]*Jc[J3+JBL];
              Ic_diffcheck[2]=-u4[0]*u4[2]*Jc[J0+JBL]+u4[2]*u4[1]*Jc[J1+JBL]+(1+u4[2]*u4[2])*Jc[J2+JBL]+u4[2]*u4[3]*Jc[J3+JBL];
              Ic_diffcheck[3]=-u4[0]*u4[3]*Jc[J0+JBL]+u4[3]*u4[1]*Jc[J1+JBL]+u4[3]*u4[2]*Jc[J2+JBL]+(1+u4[3]*u4[3])*Jc[J3+JBL];
              Is_diffcheck[0]=(1-u4[0]*u4[0])*Js[J0+JBL]+u4[0]*u4[1]*Js[J1+JBL]+u4[0]*u4[2]*Js[J2+JBL]+u4[0]*u4[3]*Js[J3+JBL];
              Is_diffcheck[1]=-u4[0]*u4[1]*Js[J0+JBL]+(1+u4[1]*u4[1])*Js[J1+JBL]+u4[1]*u4[2]*Js[J2+JBL]+u4[1]*u4[3]*Js[J3+JBL];
              Is_diffcheck[2]=-u4[0]*u4[2]*Js[J0+JBL]+u4[2]*u4[1]*Js[J1+JBL]+(1+u4[2]*u4[2])*Js[J2+JBL]+u4[2]*u4[3]*Js[J3+JBL];
              Is_diffcheck[3]=-u4[0]*u4[3]*Js[J0+JBL]+u4[3]*u4[1]*Js[J1+JBL]+u4[3]*u4[2]*Js[J2+JBL]+(1+u4[3]*u4[3])*Js[J3+JBL];
              for(l=0;l<4;l++) Ib_diffusion[l]=Jb[l+JBL]/NF-rho_b*u4[l];
              for(l=0;l<4;l++) Ic_diffusion[l]=Jc[l+JBL]/NF-rho_c*u4[l];
              for(l=0;l<4;l++) Is_diffusion[l]=Js[l+JBL]/NF-rho_s*u4[l];
              for(l=0;l<4;l++) {
                 if(fabs(Ib_diffusion[l]-Ib_diffcheck[l]/NF)>1.e-10) printf("Warning, mismatching in baryon diffusion currents at i=%d, j=%d, k=%d!  %14.11e\n",i,j,k,Ib_diffusion[l]-Ib_diffcheck[l]/NF);
                 if(fabs(Ic_diffusion[l]-Ic_diffcheck[l]/NF)>1.e-10) printf("Warning, mismatching in charge diffusion currents at i=%d, j=%d, k=%d!  %14.11e\n",i,j,k,Ic_diffusion[l]-Ic_diffcheck[l]/NF);
                 if(fabs(Is_diffusion[l]-Is_diffcheck[l]/NF)>1.e-10) printf("Warning, mismatching in strange diffusion currents at i=%d, j=%d, k=%d! %14.11e\n",i,j,k,Is_diffusion[l]-Is_diffcheck[l]/NF);
              }
              fwrite(Ib_diffusion,sizeof(double),4,fde);
              fwrite(Ic_diffusion,sizeof(double),4,fde);
              fwrite(Is_diffusion,sizeof(double),4,fde);
	      for(p=0;p<np;p++)	{
		 rho=(Jp[J0+JPL]*u4[0]-Jp[J1+JPL]*u4[1]-Jp[J2+JPL]*u4[2]-Jp[J3+JPL]*u4[3])/(cell_volume*nevents);
		 eps=((u4[0]*Tp[T00+TLOC]-u4[1]*Tp[T10+TLOC]-u4[2]*Tp[T20+TLOC]-u4[3]*Tp[T30+TLOC])*u4[0]-
		 (u4[0]*Tp[T01+TLOC]-u4[1]*Tp[T11+TLOC]-u4[2]*Tp[T21+TLOC]-u4[3]*Tp[T31+TLOC])*u4[1]-
		 (u4[0]*Tp[T02+TLOC]-u4[1]*Tp[T12+TLOC]-u4[2]*Tp[T22+TLOC]-u4[3]*Tp[T32+TLOC])*u4[2]-
		 (u4[0]*Tp[T03+TLOC]-u4[1]*Tp[T13+TLOC]-u4[2]*Tp[T23+TLOC]-u4[3]*Tp[T33+TLOC])*u4[3])/(cell_volume*nevents);
                 tmp_value=(double)Pnum[PNLOC];
                 fwrite(&tmp_value,sizeof(double),1,fde);
		 fwrite(&rho,sizeof(double),1,fde);
		 fwrite(&eps,sizeof(double),1,fde); 						 							   
	      }         
	    }
	    else
	    {
              printf("Warning, undetermined Landau frame fluid four velocity at index: i: %d, j: %d, k: %d\n",i,j,k);
              fwrite(empty_arr,sizeof(double),15+np*3,fde);
	    }						
	  
      if ((argc==4) && (success==1)){ //we compute also Tmunu in the Landau frame
          glf1=u4[0]+1;
          gsl_matrix_set(lambda_mat, 0, 0, u4[0]);
          gsl_matrix_set(lambda_mat, 0, 1, -u4[1]);
          gsl_matrix_set(lambda_mat, 0, 2, -u4[2]);
          gsl_matrix_set(lambda_mat, 0, 3, -u4[3]);
          gsl_matrix_set(lambda_mat, 1, 0, -u4[1]);
          gsl_matrix_set(lambda_mat, 1, 1, 1+u4[1]*u4[1]/glf1);
          gsl_matrix_set(lambda_mat, 1, 2, u4[1]*u4[2]/glf1);
          gsl_matrix_set(lambda_mat, 1, 3, u4[1]*u4[3]/glf1);
          gsl_matrix_set(lambda_mat, 2, 0, -u4[2]);
          gsl_matrix_set(lambda_mat, 2, 1, u4[2]*u4[1]/glf1);
          gsl_matrix_set(lambda_mat, 2, 2, 1+u4[2]*u4[2]/glf1);
          gsl_matrix_set(lambda_mat, 2, 3, u4[2]*u4[3]/glf1);
          gsl_matrix_set(lambda_mat, 3, 0, -u4[3]);
          gsl_matrix_set(lambda_mat, 3, 1, u4[3]*u4[1]/glf1);
          gsl_matrix_set(lambda_mat, 3, 2, u4[3]*u4[2]/glf1);
          gsl_matrix_set(lambda_mat, 3, 3, 1+u4[3]*u4[3]/glf1);
          
          gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,lambda_mat,Tmunu,0.,tmp_mat);
          gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,tmp_mat,lambda_mat,0.,Tmunu_Landau);
          
          //now we write the results in the output file
          fwrite(&nevents,sizeof(long int),1,ftm);
	  fwrite(&time,sizeof(double),1,ftm);
          //we save the informations about the number of particles and the grid size because, for a few bytes more of disk space, we can perform a consistency check when averaging
	  fwrite(&np,sizeof(int),1,ftm);
	  fwrite(&nx,sizeof(int),1,ftm);
      fwrite(&ny,sizeof(int),1,ftm);
      fwrite(&nz,sizeof(int),1,ftm);
      fwrite(&dx,sizeof(double),1,ftm);
      fwrite(&dy,sizeof(double),1,ftm);
      fwrite(&dz,sizeof(double),1,ftm);
      fwrite(&xmin,sizeof(double),1,ftm);
      fwrite(&ymin,sizeof(double),1,ftm);
      fwrite(&zmin,sizeof(double),1,ftm);
	  fwrite(&Tp[h*nx*ny*nz*np*10],sizeof(double),nx*ny*nz*np*10,ftm);
	  fwrite(&Jp[h*nx*ny*nz*np*4],sizeof(double),nx*ny*nz*np*4,ftm);
	  fwrite(&Jb[h*nx*ny*nz*4],sizeof(double),nx*ny*nz*4,ftm);
	  fwrite(&Jc[h*nx*ny*nz*4],sizeof(double),nx*ny*nz*4,ftm);
	  fwrite(&Js[h*nx*ny*nz*4],sizeof(double),nx*ny*nz*4,ftm);
	  fwrite(&Jt[h*nx*ny*nz*4],sizeof(double),nx*ny*nz*4,ftm);
	  fwrite(&Pnum[h*nx*ny*nz*np],sizeof(long int),nx*ny*nz*np,ftm);
	  
      }    
     }    
   }				
  }

fclose(fde);
printf("Densities data saved in file %s.\n",argv[2]);
gsl_eigen_nonsymmv_free(ws);
gsl_matrix_free(Tmunu);
gsl_matrix_complex_free(result);
gsl_vector_complex_free(eval);
if(argc==4){
gsl_matrix_free(lambda_mat);
gsl_matrix_free(tmp_mat);
gsl_matrix_free(Tmunu_Landau);
fclose(ftm);
printf("Tensor and current data saved in file %s.\n",argv[3]);
}
	   
free(Tp);
free(Jb);
free(Jp);
free(Jc);
free(Js);
free(Jt);
free(Pnum);
return 0;
}
