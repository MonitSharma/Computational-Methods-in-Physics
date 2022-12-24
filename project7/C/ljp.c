/*
   .
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Writes the coordinates in XYZ format to the output stream fp.
   The integer "z" is the atomic number of the particles, required
   for the XYZ format. The array ix contains the number of x-dir
   periodic boundary crossings a particle has performed; thus,
   the "unfolded" coordinate is rx[i]+ix[i]*L. */
void xyz_out (FILE * fp,
              double * rx, double * ry, double * rz,
              double * vx, double * vy, double * vz,
              int * ix, int * iy, int * iz, double L,
              int N, int z, int put_vel, int unfold) {
  int i;

  fprintf(fp,"%i %i\n\n",N,put_vel);
  for (i=0;i<N;i++) {
    fprintf(fp,"%i %.8lf %.8lf %.8lf ",z,
            rx[i]+(unfold?(ix[i]*L):0.0),
            ry[i]+(unfold?(iy[i]*L):0.0),
            rz[i]+(unfold?(iz[i]*L):0.0));
    if (put_vel)
      fprintf(fp,"%.8lf %.8lf %.8lf",vx[i],vy[i],vz[i]);
    fprintf(fp,"\n");
  }
}

int xyz_in (FILE * fp, double * rx, double * ry, double * rz,
             double * vx, double * vy, double * vz,
             int * N) {
  int i;
  int has_vel, dum;
  fscanf(fp,"%i %i\n\n",N,&has_vel);
  for (i=0;i<(*N);i++) {
    fscanf(fp,"%i %lf %lf %lf ",&dum,&rx[i],&ry[i],&rz[i]);
    if (has_vel) { // read velocities
      fscanf(fp,"%lf %lf %lf",&vx[i],&vy[i],&vz[i]);
    }
  }
  return has_vel;
}

/** \brief save current positions
 * \param[in] rx array of x-coordinate of the position of all particles
 * \param[in] ry array of y-coordinate of the position of all particles
 * \param[in] rz array of z-coordinate of the position of all particles
 * \param[in] N number of particles
 * \param[out] rx0 array of x-coordinate of the position of all particles
 * \param[out] ry0 array of y-coordinate of the position of all particles
 * \param[out] rz0 array of z-coordinate of the position of all particles
 */
void save_positions(double *rx,double *ry,double *rz,int N,
                    double *rx0,double *ry0,double *rz0) {
  int i;
  for (i=0;i<N;i++) {
    rx0[i]=rx[i];
    ry0[i]=ry[i];
    rz0[i]=rz[i];
  }
}

/** \brief calculates the max distance a particle has moved
 * \param[in] rx array of x-coordinate of the position of all particles
 * \param[in] ry array of y-coordinate of the position of all particles
 * \param[in] rz array of z-coordinate of the position of all particles
 * \param[in] N number of particles
 * \param[in] rx0 array of x-coordinate of the position of all particles
 * \param[in] ry0 array of y-coordinate of the position of all particles
 * \param[in] rz0 array of z-coordinate of the position of all particles
 * \return square of the max distance moved
 */
double max_moved_distance2(double *rx,double *ry,double *rz, int N,
                        double *rx0,double *ry0,double *rz0) {
  double max2=0;
  int i;
  for(i=0;i<N;i++) {
    double dx2= (rx[i]-rx0[i])*(rx[i]-rx0[i]) \
               +(ry[i]-ry0[i])*(ry[i]-ry0[i]) \
               +(rz[i]-rz0[i])*(rz[i]-rz0[i]);
    if (dx2>max2) {
      max2=dx2;
    }
  }
  return max2;
}

/** \brief calculates periodic distance between two particles
 * \param[in] i index of particle i
 * \param[in] j index of particle j
 * \param[in] rx array of x-coordinate of the position of all particles
 * \param[in] ry array of y-coordinate of the position of all particles
 * \param[in] rz array of z-coordinate of the position of all particles
 * \param[out] dx difference in x-direction
 * \param[out] dy difference in y-direction
 * \param[out] dz difference in z-direction
 * \param[in] L box length
 * \return square shorted distance of all periodic images
 */
double per_dist2(int i, int j,
    double *rx, double *ry, double *rz,
    double *dx,double *dy,double *dz,double L) {
  double  r2,hL=L/2.0;

  /* Periodic boundary conditions: Apply the minimum image
     convention; note that this is *not* used to truncate the
     potential as long as there an explicit cutoff. */

  *dx  = (rx[i]-rx[j]);
  *dy  = (ry[i]-ry[j]);
  *dz  = (rz[i]-rz[j]);

  if (*dx>hL)       *dx-=L;
  else if (*dx<-hL) *dx+=L;
  if (*dy>hL)       *dy-=L;
  else if (*dy<-hL) *dy+=L;
  if (*dz>hL)       *dz-=L;
  else if (*dz<-hL) *dz+=L;
  r2 = *dx* *dx + *dy* *dy + *dz* *dz;
  return r2;
}

/** \brief calculates Lennard-Jones interaction between two particles
 * \param[in] i index of particle i
 * \param[in] j index of particle j
 * \param[in] dx array of x-coordinate of the distance vetor
 * \param[in] dy array of y-coordinate of the distance vetor
 * \param[in] dz array of z-coordinate of the distance vetor
 * \param[in] r2 square of the distance vector
 * \param[in,out] fx array of x-coordinate of the force of all particles
 * \param[in,out] fy array of y-coordinate of the force of all particles
 * \param[in,out] fz array of z-coordinate of the force of all particles
 * \param[in] rc2 square of the cutoff distance
 * \param[in] ecut energy varlue at the cutoff
 * \param[in,out] vir variable to add the virial contribute to
 * \param[in,out] e variable to add the energy contribute to
 */
void calc_lj ( int i, int j, double dx, double dy, double dz, double r2,
               double * fx, double * fy, double * fz,
               double rc2, double ecut, double * vir, double *e) {
  double r6i, f;

  if (r2<rc2) {
    r6i   = 1.0/(r2*r2*r2);
    *e   += 4*(r6i*r6i - r6i) - ecut;
    f     = 48*(r6i*r6i-0.5*r6i);
    fx[i] += dx*f/r2;
    fx[j] -= dx*f/r2;
    fy[i] += dy*f/r2;
    fy[j] -= dy*f/r2;
    fz[i] += dz*f/r2;
    fz[j] -= dz*f/r2;
    *vir += f;
  }
}

/** \brief contains all information about the dpd thermostat */
typedef struct{
  double cutoff2; ///< square of the cutoff
  double gamma; ///< friction parameter
  double prefactor; ///< dpd constant = sqrt(24.0*T*gamma/dt)
  double *vx; ///<array of x-coordinate of the velocity of all particles
  double *vy; ///<array of y-coordinate of the velocity of all particles
  double *vz; ///<array of z-coordinate of the velocity of all particles
} dpd_t;

/** \brief calculates dpd force between two particles
 * \param[in] i index of particle i
 * \param[in] j index of particle j
 * \param[in] dx array of x-coordinate of the distance vetor
 * \param[in] dy array of y-coordinate of the distance vetor
 * \param[in] dz array of z-coordinate of the distance vetor
 * \param[in] r2 square of the distance vector
 * \param[in,out] fx array of x-coordinate of the force of all particles
 * \param[in,out] fy array of y-coordinate of the force of all particles
 * \param[in,out] fz array of z-coordinate of the force of all particles
 * \param[in] dpd structure with dpd information
 */
void calc_dpd ( int i, int j, double dx, double dy, double dz, double r2,
               double * fx, double * fy, double * fz, dpd_t* dpd) {

  double v_dot_d, amplitude;
  double tmp_x,tmp_y,tmp_z;
  if (r2 < dpd->cutoff2) {
    v_dot_d  = ( dpd->vx[i] - dpd->vx[j] ) * dx ;
    v_dot_d += ( dpd->vy[i] - dpd->vy[j] ) * dy ;
    v_dot_d += ( dpd->vz[i] - dpd->vz[j] ) * dz ;
    amplitude = dpd->prefactor / sqrt(r2) * (drand48()-0.5) - dpd->gamma * v_dot_d / r2;
    tmp_x =  amplitude * dx;
    fx[i] += tmp_x;
    fx[j] -= tmp_x;
    tmp_y =  amplitude * dy;
    fy[i] += tmp_y;
    fy[j] -= tmp_y;
    tmp_z =  amplitude * dz;
    fz[i] += tmp_z;
    fz[j] -= tmp_z;
  }
}
/** \brief contains all information about the verlet list */
typedef struct{
  int *head; ///< point to first neighboring particles of the ith particle
  int *id_head; ///< particle id of particle in head needed in combination with linked cell
  int *neighbors; ///<list of all neighboring paricles
  int size; ///< size of neighbors array
  int update; ///< if verlet list should be updated
  double cutoff; ///< verlet list cutoff
  double cutoff2; ///< square of the verlet list cutoff
  double skin2; ///<square of the verlet list skin (list cutoff - lj cutoff)
  double *rx0; ///<array of x-coordinate of the old position of all particles
  double *ry0; ///<array of y-coordinate of the old position of all particles
  double *rz0; ///<array of z-coordinate of the old position of all particles
} vlist_t;

/** \brief contains all information about the linked cell list */
typedef struct{
  int *head; ///< point to first neighboring particles of the ith cell
  int *neighbors; ///<list of all neighboring paricles
  int update; ///< if verlet list should be updated
  int ncells; ///< number of cell in each direction
  int ncells2; ///< square of number of cell in each direction
  int ncells3; ///< total of number of cells (ncells**3)
  int left; ///< if there is a neighbor on the left
  int right; ///< if there is a neighbor on the right
  double lcell; ///< length of a cell (effective cutoff)
  double cutoff; ///< list cutoff, round up from user input
  double skin2; ///<square of the linked cell list skin (list cutoff - lj cutoff)
  double *rx0; ///<array of x-coordinate of the old position of all particles
  double *ry0; ///<array of y-coordinate of the old position of all particles
  double *rz0; ///<array of z-coordinate of the old position of all particles
} lclist_t;

/** \brief contains all information about the neighbor list */
typedef struct{
  vlist_t *vl; ///< pointer to verlet list struct
  lclist_t *lc; ///< pointer to linked cell list struct
} nblist_t;

/** \brief algorithm for computing forces and potential energy.
 * algorithm for computing forces and potential energy. The virial
 *  is also computed and returned in *vir.
 */
double total_e ( double * rx, double * ry, double * rz,
                 double * fx, double * fy, double * fz,
                 int N, double L, nblist_t *nblist, dpd_t *dpd,
                 double rc2, double ecor, double ecut, double * vir ) {
   int i,j,k,l;
   int cx,cy,cz,nbx,nby,nbz,c;
   double e = 0.0;

   double dx,dy,dz,r2;

   *vir=0.0;

   /* Zero the forces */
   for (i=0;i<N;i++) {
     fx[i]=fy[i]=fz[i]=0.0;
   }

   if (!nblist) {
     /* do N square calculation */
     for (i=0;i<(N-1);i++) {
       for (j=i+1;j<N;j++) {
         r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
         calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,rc2,ecut,vir,&e);
	 if (dpd){
           calc_dpd(i,j,dx,dy,dz,r2,fx,fy,fz,dpd);
	 }
       }
     }
     return e+N*ecor;
   }

   /* rebuild cell list */
   if (nblist->lc && nblist->lc->update) {
     for(c=0;c<nblist->lc->ncells3;c++){
       nblist->lc->head[c]=-1;
     }
     /* this is an array implementation of a linked list*/
     for(i=0;i<N;i++) {
       cx=(int)floor(rx[i]/nblist->lc->lcell);
       cy=(int)floor(ry[i]/nblist->lc->lcell);
       cz=(int)floor(rz[i]/nblist->lc->lcell);
       c=cx+cy*nblist->lc->ncells+cz*nblist->lc->ncells2;
       nblist->lc->neighbors[i] = nblist->lc->head[c];
       nblist->lc->head[c] = i;
     }
     nblist->lc->update=0;
   }

   /* no verlet list, simple update */
   if(!nblist->vl){
     /* loop of cells in 3dim */
     for(cx=0;cx<nblist->lc->ncells;cx++){
       for(cy=0;cy<nblist->lc->ncells;cy++){
         for(cz=0;cz<nblist->lc->ncells;cz++){
           c=cx+cy*nblist->lc->ncells+cz*nblist->lc->ncells2;
           /* loop over all neighbor cell of c */
           for (nbx=cx-nblist->lc->left;nbx<=cx+nblist->lc->right;nbx++){
             for (nby=cy-nblist->lc->left;nby<=cy+nblist->lc->right;nby++){
                 for (nbz=cz-nblist->lc->left;nbz<=cz+nblist->lc->right;nbz++){
                 /* calc cell number with respect to periodicity */
                 int nb=((nbx+nblist->lc->ncells)%nblist->lc->ncells)
                   +((nby+nblist->lc->ncells)%nblist->lc->ncells)*nblist->lc->ncells
                   +((nbz+nblist->lc->ncells)%nblist->lc->ncells)*nblist->lc->ncells2;
                 for (i=nblist->lc->head[c];i!=-1;i = nblist->lc->neighbors[i]) {
                   for (j=nblist->lc->head[nb];j != -1;j = nblist->lc->neighbors[j]) {
                     if ( i < j ) {
                       r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
                       calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,rc2,ecut,vir,&e);
	               if (dpd){
                         calc_dpd(i,j,dx,dy,dz,r2,fx,fy,fz,dpd);
	               }
                     }
                   }
                 }
               }
             }
           }
         }
       }
     }
     return e+N*ecor;
   }

   /* verlet list update from cell lists*/
   if (nblist->lc && nblist->vl->update) {
     k=0;
     l=0;
     /* loop of cells in 3dim */
     for(cx=0;cx<nblist->lc->ncells;cx++){
       for(cy=0;cy<nblist->lc->ncells;cy++){
         for(cz=0;cz<nblist->lc->ncells;cz++){
           c=cx+cy*nblist->lc->ncells+cz*nblist->lc->ncells2;
           /* note different loop order than above*/
           for (i=nblist->lc->head[c];i!=-1;i=nblist->lc->neighbors[i]) {
             nblist->vl->head[l]=k;
             nblist->vl->id_head[l++]=i;
             /* loop over all neighbor cell of c */
             for (nbx=cx-nblist->lc->left;nbx<=cx+nblist->lc->right;nbx++){
               for (nby=cy-nblist->lc->left;nby<=cy+nblist->lc->right;nby++){
                   for (nbz=cz-nblist->lc->left;nbz<=cz+nblist->lc->right;nbz++){
                   /* calc cell number with respect to periodicity */
                   int nb=((nbx+nblist->lc->ncells)%nblist->lc->ncells)
                     +((nby+nblist->lc->ncells)%nblist->lc->ncells)*nblist->lc->ncells
                     +((nbz+nblist->lc->ncells)%nblist->lc->ncells)*nblist->lc->ncells2;
                   for (j=nblist->lc->head[nb];j != -1;j = nblist->lc->neighbors[j]) {
                     if ( i < j ) {
                       r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
                       if (r2 < nblist->vl->cutoff2) {
                         if(k == nblist->vl->size ) {
                           /* double the size */
                           nblist->vl->size*=2;
                           nblist->vl->neighbors=(int *)realloc(nblist->vl->neighbors,nblist->vl->size*sizeof(int));
                         }
                         nblist->vl->neighbors[k++]=j;
                         calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,rc2,ecut,vir,&e);
	                 if (dpd){
                           calc_dpd(i,j,dx,dy,dz,r2,fx,fy,fz,dpd);
	                 }
                       }
                     }
                   }
                 }
               }
             }
           }
         }
       }
     }
     nblist->vl->head[l]=k;
     nblist->vl->update=0;
     return e+N*ecor;
   }

   /* verlet list update without cell list */
   if (nblist->vl->update) {
     /* update neighbor list and calc lj */
     k=0;
     for (i=0;i<(N-1);i++) {
       nblist->vl->head[i]=k;
       for (j=i+1;j<N;j++) {
         r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
         if (r2 < nblist->vl->cutoff2) {
           if(k == nblist->vl->size ) {
             /* double the size */
             nblist->vl->size*=2;
             nblist->vl->neighbors=(int *)realloc(nblist->vl->neighbors,nblist->vl->size*sizeof(int));
           }
           nblist->vl->neighbors[k++]=j;
           calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,rc2,ecut,vir,&e);
	   if (dpd){
             calc_dpd(i,j,dx,dy,dz,r2,fx,fy,fz,dpd);
	   }
         }
       }
     }
     nblist->vl->head[N-1]=k;
     nblist->vl->update=0;
     return e+N*ecor;
   }

   /* use existing neighbor list with linked id mapping */
   if (nblist->vl->id_head){
     for (l=0;l<N;l++) {
       i=nblist->vl->id_head[l];
       for (k=nblist->vl->head[l];k<nblist->vl->head[l+1];k++) {
         j=nblist->vl->neighbors[k];
         r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
         calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,rc2,ecut,vir,&e);
	 if (dpd){
           calc_dpd(i,j,dx,dy,dz,r2,fx,fy,fz,dpd);
	 }
       }
     }
   } else {
     /* use existing neighbor list verlit list only */
     for (i=0;i<(N-1);i++) {
       for (k=nblist->vl->head[i];k<nblist->vl->head[i+1];k++) {
         j=nblist->vl->neighbors[k];
         r2=per_dist2(i,j,rx,ry,rz,&dx,&dy,&dz,L);
         calc_lj(i,j,dx,dy,dz,r2,fx,fy,fz,rc2,ecut,vir,&e);
	 if (dpd){
           calc_dpd(i,j,dx,dy,dz,r2,fx,fy,fz,dpd);
	 }
       }
     }
   }
   return e+N*ecor;
}

/** \brief generate gaussian random number
 * generate a gaussian random number using the
 * Box-Muller method
 * \param[in] mean mean value of the gaussian
 * \param[in] sigma standard deviation of the gaussian
 * \return gaussian random number
 */
double drand_gaussian(double mean, double sigma) {
  static int i=0;
  static double x2;
  /* http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform */
  double r,theta,u1,u2,x1;
  if (i==0) {
    i=1-i;
    u1=1-drand48();
    u2=1-drand48();
    r=sqrt(-2*log(u1))*sigma;
    theta=2*M_PI*u2;
    x1 = r*cos(theta) + mean;
    x2 = r*sin(theta) + mean ;
    return(x1);
  } else {
    i=1-i;
    return(x2);
  }
}

/* Initialize particle positions by assigning them
   on a cubic grid, then scaling positions
   to achieve a given box size and thereby, volume,
   and density */
void init ( double * rx, double * ry, double * rz,
            double * vx, double * vy, double * vz,
            int * ix, int * iy, int * iz,
            int * N, double L, double T0,
            double * KE, char * icf) {
  int i,iix,iiy,iiz;
  double cmvx=0.0,cmvy=0.0,cmvz=0.0;
  double T, fac;
  int n3=2;
  int vel_ok=0;

  /* If icf has a value, assume it is the name of a file containing
     the input configuration in XYZ format */
  if (icf) {
    FILE * fp = fopen(icf,"r");
    if (fp) {
      vel_ok = xyz_in(fp,rx,ry,rz,vx,vy,vz,N);
    }
    else {
      fprintf(stderr,"# error: could not read %s\n",icf);
      exit(-1);
    }
    fclose(fp);
  }
  /* Assign particles on a cubic lattice */
  else {

    /* Find the lowest perfect cube, n3, greater than or equal to the
       number of particles */
    while ((n3*n3*n3)<*N) n3++;

    iix=iiy=iiz=0;
    /* Assign particle positions */
    for (i=0;i<*N;i++) {
      rx[i] = ((double)iix+0.5)*L/n3;
      ry[i] = ((double)iiy+0.5)*L/n3;
      rz[i] = ((double)iiz+0.5)*L/n3;
      iix++;
      if (iix==n3) {
        iix=0;
        iiy++;
        if (iiy==n3) {
          iiy=0;
          iiz++;
        }
      }
    }
  }
  /* If no velocities yet assigned, randomly pick some */
  if (!vel_ok) {
    for (i=0;i<*N;i++) {
      vx[i]=drand_gaussian(0.0,1.0);
      vy[i]=drand_gaussian(0.0,1.0);
      vx[i]=drand_gaussian(0.0,1.0);
    }
  }
  /* Take away any center-of-mass drift; compute initial KE */
  for (i=0;i<*N;i++) {
    cmvx+=vx[i];
    cmvy+=vy[i];
    cmvz+=vz[i];
  }
  (*KE)=0;
  for (i=0;i<*N;i++) {
    vx[i]-=cmvx/(*N);
    vy[i]-=cmvy/(*N);
    vz[i]-=cmvz/(*N);
    (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  (*KE)*=0.5;
  /* if T0 is specified, scale velocities */
  if (T0>0.0) {
    T=(*KE)/(*N)*2./3.;
    fac=sqrt(T0/T);
    (*KE)=0;
    for (i=0;i<*N;i++) {
      vx[i]*=fac;
      vy[i]*=fac;
      vz[i]*=fac;
      (*KE)+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    (*KE)*=0.5;
  }
  /* Initialize periodic boundary crossing counter arrays */
  for (i=0;i<*N;i++) {
    ix[i]=0;
    iy[i]=0;
    iz[i]=0;
  }
}

int main ( int argc, char * argv[] ) {

  double * rx, * ry, * rz;
  double * vx, * vy, * vz;
  double * fx, * fy, * fz;
  int * ix, * iy, * iz;
  int N=216;
  double L=0.0;
  double rho=pow(1./1.05,3), rc2 = 1.12246, vir, vir_old, V;
  double PE, KE, TE, ecor, ecut, T0=1.0, TE0, T=1.0;
  double rr3,dt=0.001, dt2;
  int i,s;
  int nSteps = 10, fSamp=100;
  int use_e_corr=0;
  int unfold = 0;
  double langevin_gamma=0;
  double dpd_gamma=0;

  char fn[20];
  FILE * out;
  char * wrt_code_str = "w";
  char * init_cfg_file = NULL;
  double rvl=0;
  double rlc=0;

  unsigned long int Seed = 23410981;

  /* Here we parse the command line arguments;  If
   you add an option, document it in the usage() function! */
  for (i=1;i<argc;i++) {
    if (!strcmp(argv[i],"-N")) N=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dt")) dt=atof(argv[++i]);
    else if (!strcmp(argv[i],"-rc")) rc2=atof(argv[++i]);
    else if (!strcmp(argv[i],"-ns")) nSteps = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-T0")) T0=atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i],"-sd")) langevin_gamma=atof(argv[++i]);
    else if (!strcmp(argv[i],"-dpd")) dpd_gamma=atof(argv[++i]);
    else if (!strcmp(argv[i],"-fs")) fSamp=atoi(argv[++i]);
    else if (!strcmp(argv[i],"-sf")) wrt_code_str = argv[++i];
    else if (!strcmp(argv[i],"-icf")) init_cfg_file = argv[++i];
    else if (!strcmp(argv[i],"-ecorr")) use_e_corr = 1;
    else if (!strcmp(argv[i],"-seed")) Seed = (unsigned long)atoi(argv[++i]);
    else if (!strcmp(argv[i],"-rvl")) rvl = atof(argv[++i]);
    else if (!strcmp(argv[i],"-rlc")) rlc = atof(argv[++i]);
    else if (!strcmp(argv[i],"-uf")) unfold = 1;
    else if ((!strcmp(argv[i],"-h"))||(!strcmp(argv[i],"--help"))) {
      fprintf(stdout,"mdlj usage:\n");
      fprintf(stdout,"mdlj [options]\n\n");
      fprintf(stdout,"Options:\n");
      fprintf(stdout,"\t -N [integer]\t\tNumber of particles (default %i)\n",N);
      fprintf(stdout,"\t -rho [real]\t\tNumber density (default %f)\n",rho);
      fprintf(stdout,"\t -dt [real]\t\tTime step (default %f)\n",dt);
      fprintf(stdout,"\t -rc [real]\t\tCutoff radius (default %f)\n",rc2);
      fprintf(stdout,"\t -ns [real]\t\tNumber of integration steps (default %i)\n",nSteps);
      fprintf(stdout,"\t -T0 [real]\t\tInitial temperature (default %f)\n",T0);
      fprintf(stdout,"\t -sd [real]\t\tUse Langevin thermostat with frition constant 0=off (default %f)\n",langevin_gamma);
      fprintf(stdout,"\t -dpd [real]\t\tUse DPD thermostat with frition constant 0=off (default %f)\n",dpd_gamma);
      fprintf(stdout,"\t -T [real]\t\tTemperature of the Langevin thermostat(default %f)\n",T);
      fprintf(stdout,"\t -fs [integer]\t\tSample frequency (default %i)\n",fSamp);
      fprintf(stdout,"\t -sf [a|w]\t\tAppend or write config output file (default %s)\n",wrt_code_str);
      fprintf(stdout,"\t -icf [string]\t\tInitial configuration file (default %s)\n",init_cfg_file);
      fprintf(stdout,"\t -seed [integer]\tRandom number generator seed (default %li)\n",Seed);
      fprintf(stdout,"\t -rvl [real]\t\tUse neighbor lists with this cutoff (default %f)\n",rvl);
      fprintf(stdout,"\t -rlc [read]  \t\tUse linked cells with this size (default %f)\n",rlc);
      fprintf(stdout,"\t -uf          \t\tPrint unfolded coordinates in output files (default %i)\n",unfold);
      fprintf(stdout,"\t -h           \t\tPrint this info\n");
      exit(0);
    }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
              argv[i]);
      exit(-1);
    }
  }

  /* Compute the side-length */
  L = pow((V=N/rho),0.3333333);

  /* Compute the tail-corrections; assumes sigma and epsilon are both 1 */
  rr3 = 1.0/(rc2*rc2*rc2);
  ecor = use_e_corr?8*M_PI*rho*(rr3*rr3*rr3/9.0-rr3/3.0):0.0;
  ecut = 4*(rr3*rr3*rr3*rr3-rr3*rr3);

  /* Compute the *squared* cutoff, reusing the variable rc2 */
  rc2*=rc2;

  /* compute the squared time step */
  dt2=dt*dt;

  /* Output some initial information */
  fprintf(stdout,"# NVE MD Simulation of a Lennard-Jones fluid\n");
  fprintf(stdout,"# L = %.5lf; rho = %.5lf; N = %i; rc = %.5lf\n",
          L,rho,N,sqrt(rc2));
  fprintf(stdout,"# nSteps %i, seed %li, dt %.5lf\n", nSteps,Seed,dt);

  nblist_t *nblist=NULL;

  if (( rvl > 0 ) || ( rlc > 0 )){
    nblist=(nblist_t *)malloc(sizeof(nblist_t));
  }

  /* verlet list stuff */
  if (rvl > 0) {
    nblist->vl=(vlist_t *)malloc(sizeof(vlist_t));
    nblist->vl->cutoff = rvl;
    nblist->vl->cutoff2 = rvl*rvl;
    if (nblist->vl->cutoff2 < rc2) {
      fprintf(stderr,"-rvl must be bigger than -rc\n");
      exit(1);
    }
    fprintf(stdout,"# using Verlet lists with cutoff %.5f\n",nblist->vl->cutoff);
    /* assume 25 neighbors and realloc later if needed */
    nblist->vl->skin2 = pow(rvl-sqrt(rc2),2);
    nblist->vl->size= 25*N;
    nblist->vl->neighbors = (int *)malloc(nblist->vl->size*sizeof(int));
    nblist->vl->head = (int *)malloc(N*sizeof(int));
    /* always update nblist one time */
    nblist->vl->update=1;
    /* old position saver*/
    nblist->vl->rx0 = (double*)malloc(N*sizeof(double));
    nblist->vl->ry0 = (double*)malloc(N*sizeof(double));
    nblist->vl->rz0 = (double*)malloc(N*sizeof(double));
  }

  if ( rlc > 0 ) {
    nblist->lc=(lclist_t *)malloc(sizeof(lclist_t));
    nblist->lc->cutoff = rlc;
    nblist->lc->ncells=(int)floor(L/rlc);
    nblist->lc->ncells2=nblist->lc->ncells * nblist->lc->ncells;
    nblist->lc->ncells3=nblist->lc->ncells * nblist->lc->ncells2;
    nblist->lc->left=nblist->lc->right=1;
    if (nblist->lc->ncells==1){
      nblist->lc->left=nblist->lc->right=0;
    }
    else if (nblist->lc->ncells==2){
      nblist->lc->left=0;
    }
    nblist->lc->lcell = L/nblist->lc->ncells;
    nblist->lc->skin2 = pow(nblist->lc->lcell-sqrt(rc2),2);
    if (rlc*rlc < rc2) {
      fprintf(stderr,"-rlc must be bigger than -rc\n");
      exit(1);
    }
    fprintf(stdout,"# using Linked cell lists with grid %i x %i x %i (cufoff %f)\n",nblist->lc->ncells,nblist->lc->ncells,nblist->lc->ncells,nblist->lc->lcell);
    nblist->lc->neighbors = (int *)malloc(N*sizeof(int));
    nblist->lc->head = (int *)malloc(nblist->lc->ncells3*sizeof(int));
    /* always update nblist one time */
    nblist->lc->update=1;
    /* old position saver*/
    nblist->lc->rx0 = (double*)malloc(N*sizeof(double));
    nblist->lc->ry0 = (double*)malloc(N*sizeof(double));
    nblist->lc->rz0 = (double*)malloc(N*sizeof(double));
    if(nblist->vl) {
      nblist->vl->id_head = (int *)malloc(N*sizeof(int));
      /* we need to store an end head for the Nth particle, too */
      nblist->vl->head = (int *)realloc(nblist->vl->head,(N+1)*sizeof(int));
    }
  }

  if(!nblist){
    fprintf(stdout,"# no neighbor lists\n");
  }

  double langevin_const = 0;
  if (langevin_gamma > 0) {
    fprintf(stdout,"# using Langevin thermostat for T=%f with gamma=%f\n",T,langevin_gamma);
    langevin_const = sqrt(24.0*T*langevin_gamma/dt);
    if (dpd_gamma > 0 ) {
      fprintf(stderr,"# error: dpd and langevin thermostat cannot be enable at the same time\n");
      exit(-1);
    }
  }

  /* Seed the random number generator */
  srand48(Seed);

  /* Allocate the position arrays */
  rx = (double*)malloc(N*sizeof(double));
  ry = (double*)malloc(N*sizeof(double));
  rz = (double*)malloc(N*sizeof(double));

  /* Allocate the boundary crossing counter arrays */
  ix = (int*)malloc(N*sizeof(int));
  iy = (int*)malloc(N*sizeof(int));
  iz = (int*)malloc(N*sizeof(int));

  /* Allocate the velocity arrays */
  vx = (double*)malloc(N*sizeof(double));
  vy = (double*)malloc(N*sizeof(double));
  vz = (double*)malloc(N*sizeof(double));

  /* Allocate the force arrays */
  fx = (double*)malloc(N*sizeof(double));
  fy = (double*)malloc(N*sizeof(double));
  fz = (double*)malloc(N*sizeof(double));

  /* dpd thermostat */
  dpd_t *dpd = NULL;
  if(dpd_gamma > 0) {
    dpd = (dpd_t *)malloc(sizeof(dpd_t));
    dpd->gamma = dpd_gamma;
    dpd->cutoff2 = rc2; /*could be different */
    dpd->prefactor = sqrt(24.0*T*dpd_gamma/dt);
    dpd->vx = vx;
    dpd->vy = vy;
    dpd->vz = vz;
  }

  /* Generate initial positions on a cubic grid,
     and measure initial energy */
  init(rx,ry,rz,vx,vy,vz,ix,iy,iz,&N,L,T0,&KE,init_cfg_file);
  sprintf(fn,"%i.xyz",0);
  out=fopen(fn,"w");
  xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,16,1,unfold);
  fclose(out);

  if (nblist){
    if (nblist->vl) {
      save_positions(rx,ry,rz,N,nblist->vl->rx0,nblist->vl->ry0,nblist->vl->rz0);
    }
    if (nblist->lc) {
      save_positions(rx,ry,rz,N,nblist->lc->rx0,nblist->lc->ry0,nblist->lc->rz0);
    }
  }

  /* DPD thermostat */
  /* in the initial md step the thermostat only acts over "half" dt */
  if (dpd) {
    dpd->gamma = dpd_gamma * sqrt(3);
    PE = total_e(rx,ry,rz,fx,fy,fz,N,L,nblist,dpd,rc2,ecor,ecut,&vir_old);
    dpd->gamma = dpd_gamma;
  } else {
    PE = total_e(rx,ry,rz,fx,fy,fz,N,L,nblist,dpd,rc2,ecor,ecut,&vir_old);
  }
  TE0=PE+KE;

  /* Langevin thermostat */
  /* in the initial md step the thermostat only acts over "half" dt */
  if (langevin_gamma > 0) {
    for (i=0;i<N;i++) {
      fx[i]+= -langevin_gamma*vx[i] + langevin_const*sqrt(3)*(drand48()-0.5);
      fy[i]+= -langevin_gamma*vy[i] + langevin_const*sqrt(3)*(drand48()-0.5);
      fz[i]+= -langevin_gamma*vz[i] + langevin_const*sqrt(3)*(drand48()-0.5);
    }
  }

  fprintf(stdout,"# step time PotE KinE TotE Tdrift T P\n");

  for (s=0;s<nSteps;s++) {

    /* First integration half-step */
    for (i=0;i<N;i++) {
      rx[i]+=vx[i]*dt+0.5*dt2*fx[i];
      ry[i]+=vy[i]*dt+0.5*dt2*fy[i];
      rz[i]+=vz[i]*dt+0.5*dt2*fz[i];
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      /* Apply periodic boundary conditions */
      if (rx[i]<0.0) { rx[i]+=L; ix[i]--; }
      if (rx[i]>L)   { rx[i]-=L; ix[i]++; }
      if (ry[i]<0.0) { ry[i]+=L; iy[i]--; }
      if (ry[i]>L)   { ry[i]-=L; iy[i]++; }
      if (rz[i]<0.0) { rz[i]+=L; iz[i]--; }
      if (rz[i]>L)   { rz[i]-=L; iz[i]++; }
    }

    if (nblist) {
      /* auto nblist update only */
      /*2* dx > skin */
      if ((nblist->lc)&&(nblist->lc->update==0)&&
          (4*max_moved_distance2(rx,ry,rz,N,nblist->lc->rx0,nblist->lc->ry0,nblist->lc->rz0) > nblist->lc->skin2)) {
        nblist->lc->update=1;
        save_positions(rx,ry,rz,N,nblist->lc->rx0,nblist->lc->ry0,nblist->lc->rz0);
        /* trigger verlet list update too */
        if (nblist->vl){
          nblist->vl->update=1;
          save_positions(rx,ry,rz,N,nblist->vl->rx0,nblist->vl->ry0,nblist->vl->rz0);
        }
      }
      if ((nblist->vl)&&(nblist->vl->update==0)&&
          (4*max_moved_distance2(rx,ry,rz,N,nblist->vl->rx0,nblist->vl->ry0,nblist->vl->rz0) > nblist->vl->skin2)) {
        nblist->vl->update=1;
        save_positions(rx,ry,rz,N,nblist->vl->rx0,nblist->vl->ry0,nblist->vl->rz0);
      }
    }

    /* Calculate forces */
    PE = total_e(rx,ry,rz,fx,fy,fz,N,L,nblist,dpd,rc2,ecor,ecut,&vir);

    /* Langevin thermostat */
    if (langevin_gamma > 0) {
      for (i=0;i<N;i++) {
        fx[i]+= -langevin_gamma*vx[i] + langevin_const*(drand48()-0.5);
        fy[i]+= -langevin_gamma*vy[i] + langevin_const*(drand48()-0.5);
        fz[i]+= -langevin_gamma*vz[i] + langevin_const*(drand48()-0.5);
      }
    }

    /* Second integration half-step */
    KE = 0.0;
    for (i=0;i<N;i++) {
      vx[i]+=0.5*dt*fx[i];
      vy[i]+=0.5*dt*fy[i];
      vz[i]+=0.5*dt*fz[i];
      KE+=vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
    }
    KE*=0.5;
    TE=PE+KE;
    fprintf(stdout,"%i %.5lf %.5lf %.5lf %.5lf %.5le %.5lf %.5lf\n",
            s,s*dt,PE,KE,TE,(TE-TE0)/TE0,KE*2/3./N,rho*KE*2./3./N+vir/3.0/V);
    if (!(s%fSamp)) {
      sprintf(fn,"%i.xyz",!strcmp(wrt_code_str,"a")?0:s);
      out=fopen(fn,wrt_code_str);
      xyz_out(out,rx,ry,rz,vx,vy,vz,ix,iy,iz,L,N,16,1,unfold);
      fclose(out);
    }
  }
  if(nblist) {
    if(nblist->lc){
      free(nblist->lc->neighbors);
      free(nblist->lc->head);
      free(nblist->lc->rx0);
      free(nblist->lc->ry0);
      free(nblist->lc->rz0);
      free(nblist->lc);
      if(nblist->vl){
        free(nblist->vl->id_head);
      }
    }
    if(nblist->vl){
      free(nblist->vl->neighbors);
      free(nblist->vl->head);
      free(nblist->vl->rx0);
      free(nblist->vl->ry0);
      free(nblist->vl->rz0);
      free(nblist->vl);
    }
    free(nblist);
  }
  free(rx);
  free(ry);
  free(rz);
  free(vx);
  free(vy);
  free(vz);
  free(ix);
  free(iy);
  free(iz);
  free(fx);
  free(fy);
  free(fz);
  return(0);
}