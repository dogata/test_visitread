/* testing code to read in data from SILO files */

#include <silo.h>
#include <string.h>
#include <stdlib.h>
#include <complex.h>

/* declarations */
int fexists(const char *fname);
void getKdim(const char *fname, int *nkx, int *nky);
void getRdim(const char *fname, int *nx, int *ny);
void getFields(const char *fname, float **vort, float **vx, float **vy);
void writeFile(const char *fname, const int nx, const int ny, const float *field);
void makeGPIframe(const int nx, const int ny, const int nxGPI, const int nyGPI,
		  const int nxStart, const int nyStart,
		  float *fieldIn, float *fieldGPI);

int main(int argc, char *argv[])
{
  int n;
  int nx,ny,nkx,nky;
  int nxGPI = 80, nyGPI = 64; // GPI frame pixel dimensions - swapped dimension due to row-major order
  int fStart,fEnd,fNum,nFiles;
  float *vort, *vx, *vy;
  float *vortGPI, *vxGPI, *vyGPI;
  char   *fname; // filename


  printf("\nRunning script to extract from existing SILO files...\n");

  fStart = atoi(argv[1]); fEnd = atoi(argv[2]);

  printf("fStart = %d, fEnd = %d\n\n",atoi(argv[1]),atoi(argv[2]));

  //  fStart = 1734; fEnd = 1735;
  nFiles = fEnd - fStart + 1; // number of files

  // create directories
  system("mkdir vort");
  system("mkdir vx");
  system("mkdir vy");

  for (n=0;n<nFiles;n++)
    {
      fNum = fStart + n;
      asprintf(&fname,"beta.%04d.silo",fNum); // allocating sprintf

      // check for file database
      printf("fname = %s\n",fname);
      if (!fexists(fname)) { printf("\nFile does not exist! Process terminated.\n"); exit(0); }
      //file = DBOpen(fname,DB_HDF5,DB_READ);
  //if ( !file ) { printf("Error reading the file!\n"); exit(0); }

      if (n == 0)
	{
	  getRdim(fname,&nx,&ny); // get real-space dimensions
	  getKdim(fname,&nkx,&nky); // get k-space dimensions

	  vort = (float *)malloc(sizeof(float)*nx*ny); // initialize the data array
	  vx = (float *)malloc(sizeof(float)*nx*ny);
	  vy = (float *)malloc(sizeof(float)*nx*ny);
	}
  
      getFields(fname,&vort,&vx,&vy);

      //printf("vort[10] = %f\n",*(vort+10));

      free(fname);

      // extract pixels corresponding to GPI frame
      vortGPI = (float *)malloc(sizeof(float)*nxGPI*nyGPI);
      vxGPI = (float *)malloc(sizeof(float)*nxGPI*nyGPI);
      vyGPI = (float *)malloc(sizeof(float)*nxGPI*nyGPI);

      makeGPIframe(nx,ny,nxGPI,nyGPI,nx/2,ny/2,vort,vortGPI);
      makeGPIframe(nx,ny,nxGPI,nyGPI,nx/2,ny/2,vx,vxGPI);
      makeGPIframe(nx,ny,nxGPI,nyGPI,nx/2,ny/2,vy,vyGPI);

      // write to files
      asprintf(&fname,"vort/vort.%04d.txt",fNum);
      writeFile(fname,nxGPI,nyGPI,vortGPI);
      free(fname);

      asprintf(&fname,"vx/vx.%04d.txt",fNum);
      writeFile(fname,nxGPI,nyGPI,vxGPI);
      free(fname);

      asprintf(&fname,"vy/vy.%04d.txt",fNum);
      writeFile(fname,nxGPI,nyGPI,vyGPI);
      free(fname);

    }


  // clear memory
  free(vort);
  free(vx);
  free(vy);

  printf("\nProcess ended.\n");
  
  return 0;
}

/* Check for file */
int fexists(const char *fname)
{
  FILE *file;
  if (file = fopen(fname,"r"))
  {
    fclose(file);
    return 1;
  }
  return 0;
}


/* get the vorticity, vx, and vy fields
   -- since the pointeresses are being altered, need to pass in type (float**)
*/
void getFields(const char *fname, float **vort, float **vx, float **vy)
{
  int n,ntot;
  DBfile *fileDB = NULL; // file pointer
  DBquadvar *vortDB = NULL; // variable pointer
  DBquadvar *vxDB = NULL, *vyDB = NULL;

  fileDB = DBOpen(fname,DB_HDF5,DB_READ);   // open database
  
  // read variable
  vortDB = DBGetQuadvar(fileDB,"vorticity");
  //printf("phi = %s\n",(*phiDB).name);
  vxDB = DBGetQuadvar(fileDB,"velocityX");
  vyDB = DBGetQuadvar(fileDB,"velocityY");

  ntot = (*vortDB).nels;
  //printf("datatype = %d\n",(*vortDB).datatype);
  //printf("nels = %d\n",(*vortDB).nels);
  //printf("nvals = %d\n",(*vortDB).nvals);

  *vort = ((float*)*((*vortDB).vals)); // syntax to dereference double void pointer in a structure!
  *vx = ((float*)*((*vxDB).vals));
  *vy = ((float*)*((*vyDB).vals));

  //if ( (*vortDB).major_order == 1 ) { printf("row-major order\n"); }
  //else { printf("column-major order\n"); }
  //printf("major_order = %d\n",(*vortDB).major_order);

  // printf("((float*)*((*vortDB).vals))[10] = %f\n",n,((float*)*((*vortDB).vals))[10]); 
  // printf("*vort[10] = %f\n",*(vort+10));

  DBClose(fileDB);
  //  return 0;
}


/* get dimension in k-space */
void getKdim(const char *fname, int *nkx, int *nky)
{

  DBfile *fileDB = NULL;
  DBquadmesh *kmesh = NULL;

  fileDB = DBOpen(fname,DB_HDF5,DB_READ); // open database
  kmesh = DBGetQuadmesh(fileDB,"kmesh"); // get mesh data
    
  *nkx = (*kmesh).dims[0];
  *nky = (*kmesh).dims[1];

  DBClose(fileDB);
}

/* get dimension in real-space */
void getRdim(const char *fname, int *nx, int *ny)
{

  DBfile *fileDB = NULL;
  DBquadmesh *smesh = NULL;

  fileDB = DBOpen(fname,DB_HDF5,DB_READ); // open database
  smesh = DBGetQuadmesh(fileDB,"spacemesh"); // get mesh data
    
  *nx = (*smesh).dims[0];
  *ny = (*smesh).dims[1];

  DBClose(fileDB);
}


/* write the field to file 
   -- to write out the field, the pointer doesn't need to be altered so pointer is of type (float*)
*/
void writeFile(const char *fname, const int nx, const int ny, const float *field)
{
  int n,m;
  FILE *fout;

  //printf("field[10] = %f\n",field[10]);

  fout = fopen(fname,"w+");
  for (n = 0; n < nx; n++)
    {
      for (m = 0; m < ny; m++)
	{
	  fprintf(fout,"%f ",field[n*ny+m]);
	}
      fprintf(fout,"\n");
    }

}

/* extracts the GPI frame */
void makeGPIframe(const int nx, const int ny, const int nxGPI, const int nyGPI,
		  const int nxStart, const int nyStart,
		  float *fieldIn, float *fieldGPI)
{
  int n,m,p,q;

  p = 0;

  //  printf("nxStart = %d, nyStart = %d\n",nxStart,nyStart);
  //  printf("nxGPI = %d, nyGPI = %d\n",nxGPI,nyGPI);

  // require to loop through all dimensions
  for (n = 0; n < nx; n++)
    {
      if ((n >= nxStart) && (n < nxStart+nxGPI))
	{
	  q = 0;
	  for (m = 0; m < ny; m++)
	    {
	      if ((m >= nyStart) && (m < nyStart+nyGPI))
		{
		  //printf("fieldIn[%d*%d+%d] = %f\n",n,nx,m,fieldIn[n*nx+m]);
		  fieldGPI[p*nyGPI + q] = fieldIn[n*ny+m];
		  //printf("fieldGPI[%d*%d + %d] = %f\n",p,nxGPI,q,fieldGPI[p*nxGPI + q]);
		  q++;
		  //printf("q = %d\n",q);
		}
	    }
	  p++;
	  //printf("p = %d\n",p);
	}
    }
}

/* //Perform FFTW
void doFFT(const int nx, const int ny, real )
{
  //fftw_real rout;
  //fftw_complex cin;
  rfftwnd_plan planb;
  
  // create plan
  planb = rfftw2d_create_plan(nx,ny,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);

  // perform FFT
  rfftwnd_one_complex_to_real(planb,cin,rout);
  
  //destroy plan
  rfftwnd_destroy_plan(planb);
}
*/
