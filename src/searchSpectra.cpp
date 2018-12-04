#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int findMinIdx(double d[], int len) {
  int iMin = 0;
  for (int i = 1; i < len; i++){
    if (d[iMin] > d[i]) {
      iMin = i;
    }
  }
  // cout<< iMin <<endl;
  return iMin;
}

// [[Rcpp::export]]
List searchSpectra(List specExp, NumericMatrix mzRange, NumericVector mz) {
  int lenMZ = mz.length();
  int lenSpec = specExp.size();
  Rcpp::List result(lenMZ);
  
  // double specMatched[lenMZ][lenSpec][2];
  for (int iMZ = 0; iMZ < lenMZ; iMZ++) {
    NumericVector mzSpec(lenSpec), intSpec(lenSpec);
    
    for (int iSpec = 0; iSpec < lenSpec; iSpec++) {
      NumericMatrix spec = as <NumericMatrix> (specExp[iSpec]);
      double mz1 = mz(iMZ);
      // cout << mz1 << '\t';
      double mzStart = mzRange(iMZ, 0);
      double mzEnd = mzRange(iMZ, 1);
      double mzDiff[spec.nrow()];
      int nr = spec.nrow();
      
      for (int iRec = 0; iRec < nr; iRec++) {
        mzDiff[iRec] = abs(spec(iRec, 0) - mz1);
        // cout << mzDiff[iRec] << '\t';
      }
      int iMin = findMinIdx(mzDiff, nr);
      double mzClosest = spec(iMin, 0);
      if (mzClosest >= mzStart && mzClosest <= mzEnd) {
        mzSpec(iSpec) = mzClosest;
        intSpec(iSpec) = spec(iMin, 1);
      }
      else {
        mzSpec(iSpec) = mz1;
        intSpec(iSpec) = 0;
      }
      // cout << mzStart << '\t' << spec(iMin-1, 0) << '\t' << mzClosest << '\t' << spec(iMin+1, 0) << '\t' << mzEnd << endl;
      // cout << mzSpec(iSpec) << '\t' << intSpec(iSpec) << endl;
    }
    result(iMZ) = Rcpp::List::create(Rcpp::Named("mz") = mzSpec,
           Rcpp::Named("intensity") = intSpec);
  }
  return result;
}
