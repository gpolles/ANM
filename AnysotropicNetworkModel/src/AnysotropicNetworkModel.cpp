//============================================================================
// Name        : AnysotropicNetworkModel.cpp
// Author      : Guido Polles
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include <mylib/mylib.h>

using namespace std;
using namespace mylib;

int main(int argc, char** argv) {
  PdbFile pdb(argv[1]); // already only CA, no check

  vector<Vector3d> crd = pdb.coords;

  ContactMap cm(crd,7.5);

  SparseMatrix<float> H0(crd.size()*3);





	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}

void generate_matrix(vector<Vector3d>& crd, ContactMap& cm, SparseMatrix<float>& H){
  size_t N = crd.size();
  for (size_t i = 0; i < N; ++i){

    Matrix3d diagonal;

    for (size_t j : cm[i]){
      if (j<i) continue;
      Vector3d dist = crd[i]-crd[j];
      double d2=normSQ(dist);
      for (int alpha = 0; alpha < 3; ++alpha) {
        for (int beta = 0; beta < 3; ++beta) {
          double v = dist[alpha]*dist[beta]/d2;
          H.insert(i*3+alpha,j*3+beta, -v);
          H.insert(j*3+beta,i*3+alpha, -v);
          diagonal[alpha][beta]+=v;
        }
      }
    }

    // insert "diagonal" elements
    for (int alpha = 0; alpha < 3; ++alpha) {
      for (int beta = 0; beta < 3; ++beta) {
        H.insert(i*3+alpha,i*3+beta,diagonal[alpha][beta]) ;
      }
    }
  }
}

void integrate_dof(vector<size_t>& intids, ContactMap& newcm, SparseMatrix<float>& H){
  for (size_t id : intids){
    for (size_t alpha = 2; alpha >= 0; --alpha) {
      size_t k = id*3+alpha;
      double l_k = 1.0/H(k,k);
      auto& row = H.getRow(k);

      // apply corrections
      for (size_t ri = 0; ri < row.size(); ++ri) {
        for (size_t rj = ri+1; rj < row.size(); ++rj) {
          size_t m = row[ri].columnIndex;
          size_t n = row[rj].columnIndex;
          float H_km = row[ri].value;
          float H_kn = row[rj].value;
          float correction = l_k*H_km*H_kn;
          // correct only for i,j in contact
          size_t i = m/3;
          size_t j = m/3;
          if ( newcm.isNeighbor(i,j) ){
            if (!H.isDefined(m,n)){
              H.insert(m,n,correction);
              H.insert(n,m,correction);
            }else{
              H(m,n) += correction;
              H(n,m) += correction;
            }
          }
        }
      }

      // remove k-row and column
      for (auto el : row){
        size_t col = el.columnIndex;
        H.remove(col,k);
      }
      std::vector<SparseMatrixElement<float> >().swap(row);
    }
  }

  // recalculate diagonals
  for (size_t i : intids){
    Matrix3d diagonal;
    for (int alpha = 0; alpha < 3; ++alpha) {
      size_t k = i*3+alpha;
      auto& row = H.getRow(k);
      for (auto el: row){
        size_t j = el.columnIndex / 3;
        size_t beta = el.columnIndex % 3;
        if (i!=j){
          diagonal[alpha][beta] -= el.value;
        }
      }
    }
    for (size_t alpha = 0; alpha < 3; ++alpha) {
      for (size_t beta = 0; beta < 3; ++beta) {
        H(i*3+alpha,i*3+beta) = diagonal[alpha][beta];
      }
    }
  }

}
