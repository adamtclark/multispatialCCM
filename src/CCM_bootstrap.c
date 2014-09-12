#include <R.h>
#include <Rmath.h>
#include <stdio.h>
int getorder(int *neighbors, double distances[], int E, int from, int to, int i, int LibUse[]);
void getrho(double A[], double Aest[], double rho[], int from, int to, int l, int *acceptablelib, int lengthacceptablelib);
void CCM_bootstrap(double *A, double *Aest, double *B, double *rho, int *pE, int *ptau, int *plengtht, int *pLibLength, int *DesiredL, int *plengthDesL, int *acceptablelib, int *plengthacceptablelib) {
    int i, j, k, l, from, to, slide, count, lindex;
    double distsv, sumu, sumaest, sumw;
    int E = *pE;
	  int nneigh=E+1;
    int tau = *ptau;
    int lengtht= *plengtht;
    int lengthDesL= *plengthDesL;
    int LibLength= *pLibLength;
    int neighbors[E+1];
    double u[E+1], w[E+1], distances[LibLength];
    int LibUse[LibLength];
    int lengthacceptablelib = *plengthacceptablelib;
    int integerpos;
    
    GetRNGstate(); // Load seed for random number generation, R base
    
    /* Code to implement Sugihara&al 2012 CCM algorithm to determing causality */
    /* Checks to see that *A causes *B */
    /* Note that *A and *B must be the same length, and standardized to same timestep */
    /* *distances is length(*A), and *Aest are length (length(*A or *B)- tau -(E+1)) */
    /* *E, *tau, and *lengtht are all single integers */
    /* *neighbors, *u, and *w are all (E+1) long */
    
    from=(tau*(E-1)); // starting point for most vectors
    for(lindex=0; lindex<lengthDesL; lindex++) { /* try out all desired library sizes, within ((E+1) to LibLengh-tau*(E-1)) */
        l=DesiredL[lindex];
        if(l<(from+(E+1))) {
            l=(from+(E+1));
        } /* Catch values that fall outside of feasible library */
        if(l>=lengtht) {
            l=(lengtht-1);
        } /* Catch values that fall outside of feasible library */
        to=l; // Set "end" of library for each iteration

			integerpos=floor(runif(0,1)*(lengthacceptablelib));
      count=acceptablelib[integerpos]; //Random number generator - populates Count with an entry that has enough of a history for given E.
      			
			for(j=from; j<=to; j++){ //For each desired library length, pull random, "acceptable" elements from library
					LibUse[j]=count;
         	integerpos=floor(runif(0,1)*(lengthacceptablelib));
      		count=acceptablelib[integerpos]; //Random number generator
			} // LibUse is now a vector of positions, based on acceptable starting points (have sufficient lags for neighbor and prediction, don't "jump" over gaps)
                    
                    for(int ii=0; ii<lengthacceptablelib; ii++) { //Predict all points in A using information from the chosen
                    	i=acceptablelib[ii]; // Make sure we only predict points that have suitable time lag information
                        for(j=from; j<=to; j++) { // scroll across elements in minimized L (based on lengthDesL, including wrapping)
                            distances[LibUse[j]]=0;
                            for(k=0; k<E; k++) {
                                distances[LibUse[j]]=distances[LibUse[j]]+pow((B[i-tau*k]-B[LibUse[j]-tau*k]),2); //calculate distances between focul point i and all other points j on shadow manifold for all E lagged dimensions
                            }
                            distances[LibUse[j]]=sqrt(distances[LibUse[j]]);
                        }
                        //distances is now a vector of distances between focul point, and all j other points in simulated library under consideration, and is indexed by LibUse[]
                        nneigh = getorder(neighbors, distances, E, from, to, i, LibUse); //find position of (E+1) closest points to B[i] on the shadow manifold. Neighbors correspond to positions in LibUse.
                        distsv=distances[neighbors[0]]; //shortest distance
                        
                        sumaest=0.; /* find w, and weighted Aest variables */
                        if(distsv!=0) { /* check whether minimum distance is zero */
                            sumu=0.; /* find u for all neighbors, and sumu */
                            for(j=0; j<(nneigh); j++) {
                                u[j]=exp(-distances[neighbors[j]]/distsv);
                                sumu=sumu+u[j];
                            }
                            sumw=0.;
                            for(j=0; j<(nneigh); j++) {
                                w[j]=u[j]/sumu;
                                if(w[j]<0.000001) { //minimum cap on weights, taken from Hao Ye code (Sugihara CCM paper)
                                    w[j]=0.000001;
                                }
                                sumw=sumw+w[j];
                            }
                            for(j=0; j<(nneigh); j++) {
                                w[j]=w[j]/sumw;
                                sumaest=sumaest+(A[neighbors[j]])*(w[j]);
                            }
                        } else {
                            sumw=0.;
                            for(j=0; j<(nneigh); j++) {
                                if(distances[neighbors[j]]==0) {
                                    w[j]=1;
                                } else {
                                    w[j]=0.000001;
                                }
                                sumw=sumw+w[j];
                            }
                            for(j=0; j<(nneigh); j++) {
                                w[j]=w[j]/sumw;
                                sumaest=sumaest+(A[neighbors[j]])*(w[j]);
                            }
                        }
                        Aest[i]=sumaest; /* calculate Aest */
                    }
                    getrho(A, Aest, rho, from, to, l, acceptablelib, lengthacceptablelib);
            rho[l]=rho[l]; // Calcualte mean of rho values
    } // End desired libraries
    PutRNGstate(); // Free up state of R random number generator
}

int getorder(int *neighbors, double distances[], int E, int from, int to, int i, int LibUse[]) {

    //find position of (Ego+1) closest points to B[i] on the shadow manifold
    int trip, n, ii, j, k;
    int nneigh=1;
    n=0;
    if(LibUse[from]==i) { //if first element is a self-reference, increment
        n=1;
    }
	  neighbors[0]=LibUse[from+n];
	
    for(ii=(from+n); ii<=to; ii++) { // scroll across elements in L
        trip=0;
        for(j=0; j<nneigh; j++) { // move existing neighbors up in vector if a closer neighbor is found
            if((distances[LibUse[ii]]<distances[neighbors[j]])&&(LibUse[ii]!=i)) { // check whether any distance is smaller than distances[neighbors[]]
                for(k=(nneigh); k>(j); k--) { // Move all remaining neighbors up a position in the vector
                    if(k<(E+1)) {
                        neighbors[k]=neighbors[k-1];
                    }
                }
                neighbors[j]=LibUse[ii]; // Replace leading element with new closest neighbor

                if((nneigh<(E+1))&&(trip!=0)) {
                    nneigh=nneigh+1;
                }
                trip=1;
                break;
            }
        }
        if((trip==0)&&(nneigh<(E+1))&&(LibUse[ii]!=i)&&(neighbors[nneigh-1]!=LibUse[ii])) { // Add element if vector is not yet full, but distance is greater than all already recorded
            neighbors[nneigh]=LibUse[ii];
            if(nneigh<(E+1)) {
				      nneigh=nneigh+1;
            }
        }
    }
    return nneigh;
}

void getrho(double A[], double Aest[], double rho[], int from, int to, int l, int *acceptablelib, int lengthacceptablelib) {
    // Calculate Pearson correlation coefficient between A and Aest
    int j, jj, n=0;
    double xbar=0, ybar=0, rhocalc=0, xyxybar=0, xxbarsq=0, yybarsq=0;

    for(jj=0; jj<lengthacceptablelib; jj++) { // scroll across elements that were used in acceptable library, and calculate means for A
        j=acceptablelib[jj];
        xbar=xbar+A[j];
        ybar=ybar+Aest[j];
        n=n+1;
    }
    xbar=xbar/n;
    ybar=ybar/n;
    
    for(jj=0; jj<lengthacceptablelib; jj++) { // Calculate mean residuals
        j=acceptablelib[jj];
        xyxybar=xyxybar+((A[j]-xbar)*(Aest[j]-ybar));
        xxbarsq=xxbarsq+pow(A[j]-xbar,2);
        yybarsq=yybarsq+pow(Aest[j]-ybar,2);
    }
    rhocalc=xyxybar/(sqrt(xxbarsq)*sqrt(yybarsq));
    if((rhocalc>=-1)&&(rhocalc<=1)) {
      rho[l]=rho[l]+rhocalc;
    }
}

