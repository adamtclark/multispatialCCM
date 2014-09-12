#include <R.h>
#include <Rmath.h>
#include <stdio.h>
void getorder_ssr(int neighbors[], double distances[], int E, int lengthacceptablelib, int i, int acceptablelib[], int predstep, int nneigh);
void SSR_predict_boot(double *A, double *Aest, double *B, int *pE, int *ptau, int *pBlength, int *pAlength, int *ppredstep, int *prepvec, int *pmatchSugi, int *acceptablelib, int *plengthacceptablelib) {
    int i=0, j=0, k=0, from, ii=0, jj=0;
    double distsv=0, sumu, sumaest, sumw;
    int E = *pE;
    int nneigh=E+1;
    int tau = *ptau;
    int Blength= *pBlength;
    int Alength= *pAlength;
    int repvec= *prepvec;
    int neighbors[E+1];
    int predstep= *ppredstep;
    int matchSugi=*pmatchSugi;
    double u[E+1], w[E+1], distances[Blength];
    double maxdist=0;
    int lengthacceptablelib=*plengthacceptablelib;
    /* Code to implement Sugihara&al SSR algorithm */
    /* Uses information in B to predict next predstep steps of A */
    /* Note that *A and *B need not be the same size - but if they are the same vector, then repvec must = "1" */
    
    
    from=(tau*(E-1)); /* starting point in most vectors */
    for(ii=0; ii<lengthacceptablelib; ii++) {	/*  include all indices  */
    	i = acceptablelib[ii]; // Consider only elements with sufficient information for lags
        if(repvec==1) { /* A and B are same vector - use "leave one out cross validation" */
            for(jj=0; jj<(lengthacceptablelib-predstep); jj++) { /* scroll across elements in L */
                j = acceptablelib[jj]; // Consider only elements with sufficient information for lags

                distances[j]=0;

                if(matchSugi==1) {
                if(i!=j) {
                    for(k=0; k<E; k++) {
                        distances[j]=distances[j]+pow((A[i-tau*k]-B[j-tau*k]),2); /*calculate distances between points on A and B for all E lagged dimensions*/
                    }
                    distances[j]=sqrt(distances[j]);
                    if(maxdist<distances[j]) {
                      maxdist=999999999*distances[j];
                    }
                } else {
                    distances[j]=maxdist;
                }
                } else{
                    if((j>i+predstep)|(j<=(i-E))) {
                        for(k=0; k<E; k++) {
                            distances[j]=distances[j]+pow((A[i-tau*k]-B[j-tau*k]),2); /*calculate distances between points on A and B for all E lagged dimensions*/
                        }
                        distances[j]=sqrt(distances[j]);
                        if(maxdist<distances[j]) {
                          maxdist=999999999*distances[j];
                        }
                    } else {
                        distances[j]=maxdist;
                    }
                }
            }
        } else {
            for(jj=0; jj<lengthacceptablelib; jj++) { /* scroll across elements in L */
            	j = acceptablelib[jj]; // Consider only elements with sufficient information for lags
                distances[j]=0;
                for(k=0; k<E; k++) {
                    distances[j]=distances[j]+pow((A[i-tau*k]-B[j-tau*k]),2); /*calculate distances between points on A and B for all E lagged dimensions*/
                }
                distances[j]=sqrt(distances[j]);
            }
        }
	getorder_ssr(neighbors, distances, E, lengthacceptablelib, i, acceptablelib, predstep, nneigh); /*find "tme" position of (E+1) closest points to A[i] on B manifold*/
        distsv=distances[neighbors[0]];

if(nneigh==(E+1)) { // Make sure we found enough neighbors

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
                if(w[j]<0.000001) {
                    w[j]=0.000001;
                }
                sumw=sumw+w[j];
            }
            for(j=0; j<(nneigh); j++) {
                w[j]=w[j]/sumw;
                sumaest=sumaest+(B[neighbors[j]+predstep])*(w[j]);
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
} else {
  Aest[i]=0;
  } // end if(nneigh == E+1)

    }
}


void getorder_ssr(int neighbors[], double distances[], int E, int lengthacceptablelib, int i, int acceptablelib[], int predstep, int nneigh) {
    /*find "tme" position of (Ego+1) closest points to B[i] on the shadow manifold*/
    /* include output only for distancesgo >0 */
    int trip, n, ii, iii, j, k;    
    nneigh=1;
    n=0;

    if(acceptablelib[0]==i) { /* if first element is a self-reference */
        n=n+1; /* #### scroll across elements in L (including wrapping) #### */
    }
    if(n>=lengthacceptablelib) {
    	n=lengthacceptablelib-1;
    }
    neighbors[0]=acceptablelib[0+n];
    for(iii=(0+n); iii<lengthacceptablelib-predstep; iii++) { /* #### scroll across elements in L (including wrapping) #### */
     	ii=acceptablelib[iii];
        trip=0;
        for(j=0; j<nneigh; j++) {
            if((distances[ii]<distances[neighbors[j]])&(ii!=i)&(j>0)) {
                for(k=(nneigh); k>(j); k--) { // slide all elements back
                    if(k<(E+1)) {
                        neighbors[k]=neighbors[k-1];
                    }
                }
                neighbors[j]=ii;
                if(nneigh<(E+1)) {
                    nneigh=nneigh+1;
                }
                trip=1;
                break;
            }
        }
        if((trip==0)&(nneigh<(E+1))&(ii!=i)&(neighbors[nneigh-1]!=ii)) { /* Add element if vector is not yet full, but distance is greater than all already recorded */
            neighbors[nneigh]=ii;
            if(nneigh<(E+1)) {
				      nneigh=nneigh+1;
            }
        }
    }
}
