/*
  haldane.c

  Solve basic DE for Haldane model with piecewise linear pressure function

  $Revision: 1.5 $ $Date: 2022/10/18 09:58:05 $

 */

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

/* how to address a matrix */
#define MATRIX(M,I,J,NI) M[I + NI * J]
#define ARRAY(A,I,J,K,NI,NJ) A[I + NI * (J + NJ * K)]


void HaldaneCalc(int    *ntimes,
		 int    *ntissues,
		 int    *ngases, 
		 double *initial,
		 double *pressures,
		 double *durations,
		 double *gasfractions,
		 double *ratecoefs,
		 int    *progressive,
		 double *profile,
		 double *final) 
{
  int itime, jtissue, kgas, Ntimes, Ntissues, Ngases, Ntime1;
  double *state;
  double p0, p1, dp, tim, fgas, r, v;

  Ntimes   = *ntimes;
  Ntissues = *ntissues;
  Ngases   = *ngases;

  /* allocate space for current state vector */
  state = (double *) R_alloc(Ntissues * Ngases, sizeof(double));

#define STATE(J,K) MATRIX(state, J, K, Ntissues)
#define INITIAL(J,K) MATRIX(initial, J, K, Ntissues)
#define RATECOEF(J,K) MATRIX(ratecoefs, J,K, Ntissues)
#define GASFRACTION(I,K) MATRIX(gasfractions, I,K, Ntimes)
#define PROFILE(I,J,K) ARRAY(profile, I,J,K, Ntimes, Ntissues)
#define FINAL(J,K) MATRIX(final, J,K, Ntissues)

  /* initialise current state */
  for(jtissue = 0; jtissue < Ntissues; jtissue++) 
    for(kgas = 0; kgas < Ngases; kgas++)
      STATE(jtissue, kgas) = INITIAL(jtissue, kgas);

  if(*progressive == 1) {
    /* copy initial state into profile */
    for(jtissue = 0; jtissue < Ntissues; jtissue++) 
      for(kgas = 0; kgas < Ngases; kgas++)
	PROFILE(0, jtissue, kgas) = INITIAL(jtissue, kgas);
  }

  /* run algorithm */
  if(Ntimes > 1) {
    if(*progressive == 1) {
      Ntime1 = Ntimes - 1;
      for(itime = 0; itime < Ntime1; itime++) {
	tim = durations[itime];
	if(tim > 0) {
	  p0 = pressures[itime];
	  p1 = pressures[itime + 1];
	  dp = p1 - p0;
	  for(kgas = 0; kgas < Ngases; kgas++) {
	    fgas = GASFRACTION(itime, kgas);
	    for(jtissue = 0; jtissue < Ntissues; jtissue++) {
	      r = RATECOEF(jtissue, kgas);
	      v = exp(- r * tim);
	      /* copy current state into profile */
	      PROFILE((itime+1), jtissue, kgas) =
		STATE(jtissue, kgas) =
		STATE(jtissue, kgas) * v + 
		fgas * (dp + (1-v) * (p0 - dp/(r*tim)));
	    }
	  }
	}
      }
    } else {
      Ntime1 = Ntimes - 1;
      for(itime = 0; itime < Ntime1; itime++) {
	tim = durations[itime];
	if(tim > 0) {
	  p0 = pressures[itime];
	  p1 = pressures[itime + 1];
	  dp = p1 - p0;
	  for(kgas = 0; kgas < Ngases; kgas++) {
	    fgas = GASFRACTION(itime, kgas);
	    for(jtissue = 0; jtissue < Ntissues; jtissue++) {
	      r = RATECOEF(jtissue, kgas);
	      v = exp(- r * tim);
	      STATE(jtissue, kgas) =
		STATE(jtissue, kgas) * v + 
		fgas * (dp + (1-v) * (p0 - dp/(r*tim)));
	    }
	  }
	}
      }
    }
  }
  /* copy final state */
  for(kgas = 0; kgas < Ngases; kgas++) 
    for(jtissue = 0; jtissue < Ntissues; jtissue++) 
      FINAL(jtissue, kgas) = STATE(jtissue, kgas);
}
