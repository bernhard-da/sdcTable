#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void calcMinimum(int *vec, int *erg, int *nDims) {
  int i, j, anz, to, z;
  j = 0;
  anz = (int)pow(2, nDims[0]);
  for (i=0; i < nDims[0]; i++) {
    to = j + anz;
    for (z = j; z < to; z++) {
      j = j + 1;
      if (j % anz == 0) {
        erg[i] = vec[z];
      }
      else {
        if(vec[z] < erg[i]) {
          erg[i] = vec[z];
        }
      }
    }
  }
}

void calcAggregationsstufen(int *indices, int *minDims, int *nDims) {
  int i, j, anz, to, z;
  anz = (int)pow(2, nDims[0]);
  j = 0;

  for (i = 0; i < nDims[0]; i++) {
    to = j + anz;
    for (z = j; z < to; z++) {
      j = j + 1
        if(indices[j] == minDims[i]) {
          indices[z] = 2;
        }
        else {
          indices[z] = 1;
        }
    }
  }
}

void calcQuader(int *g, int *dia, int *nDims, int *out) {
  int i, j, k, each, t, z, anzGruppen;
  int tmperg;
  j = 0;

  for(i=0; i < nDims[0]; i++) {
    each = (int)pow(2, (i));
    anzGruppen = (int)pow(2,nDims[0])/each;
    int erg[anzGruppen];

    for (t=0; t < anzGruppen; t++) {
      tmperg = (t+1) % 2;
      erg[t] = tmperg;
    }

    for (k=0; k < anzGruppen; k++)  {
      for(z=1; z <= each; z++) {
        if(erg[k] != 0) {
          out[j] = g[i];
        }
        else {
          out[j] = dia[i];
        }
        j++;
      }
    }
  }
}

void calcQuaderPosition(int *vals, int *lenVals, int *valsQ, int *erg, int *nDims) {
  int i, j, indikator;
  int anz = (int)pow(2, nDims[0]);
  indikator = 1;
  for (i=0; i < anz; i++) {
    j = 0;
    while(indikator > 0) {
      if(valsQ[i] == vals[j]) {
        erg[i] = j+1;
        indikator = -1;
      }
      else {
        j++;
      }
    };
    indikator = 1;
  }
}

void extractIndicesSubtable(int *vec, int *lengthSub, int *erg, int *nDims, int *powers, int *final) {
  int summe = 0;
  int i, j;
  int lengthSubVec = lengthSub[0];
  int index;

  /* Wir berechnen die notwendigen Potenzen! */
  for (i=1; i <= nDims[0]; i++) {
    for (j=1; j <= lengthSubVec; j++) {
      index = ((i-1)*lengthSubVec)+j-1;
      if(vec[index] > erg[(i-1)]) {
        erg[(i-1)] = vec[index];
      }
    }
    erg[(i-1)] = (int)(log10(erg[(i-1)])+1);
    summe = summe + erg[(i-1)];
  }

  /* potenzen generieren */
  powers[0] = summe - 1;
  for (i=1; i < nDims[0]; i++) {
    powers[i] = powers[(i-1)] - erg[i];
  }

  /* Wir berechnen den Rest */
  summe = 0;
  int k = 0;
  for (j=1; j <= lengthSubVec; j++) {
    for (i=1; i <= nDims[0]; i++) {
      index = ((i-1)*lengthSubVec)+j-1;
      summe = summe + vec[index]* (int)pow(10, powers[i-1]);
    }
    final[k] = summe;
    k++;
    summe = 0;
  }
}

void extractIndicesAktQuader(int *vec, int *lengthSub, int *nDims, int *powers, int *final) {
  int summe = 0;
  int i, j;
  int lengthSubVec = lengthSub[0];
  int index;

  /* Wir berechnen den Rest */
  int k = 0;
  for (j=1; j <= lengthSubVec; j++) {
    for (i=1; i <= nDims[0]; i++) {
      index = ((i-1)*lengthSubVec)+j-1;
      summe = summe + vec[index]* pow(10, powers[i-1]);
    }
    final[k] = summe;
    k++;
    summe = 0;
  }
}

void normQuader(int *indices, int *nDims, int *lVec) {
  int i, j, index;
  int N     = (lVec[0] / nDims[0]);

  for (i=2; i <= N; i++) {
    for (j=0; j < nDims[0]; j++) {
      index = (i-1)*nDims[0]+j;
      if(indices[index] == indices[j]) {
        indices[index] = 0;
      }
      else {
        indices[index] = 1;
      }
    }
  }
  for (j=0; j < nDims[0]; j++) {
    indices[j] = 0;
  }
}
