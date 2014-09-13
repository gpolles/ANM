//============================================================================
// Name        : AnysotropicNetworkModel.cpp
// Author      : Guido Polles
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <mylib/mylib.h>

using namespace std;
using namespace mylib;

template<class T>
void generate_matrix(vector<Vector3d>& crd, ContactMap& cm, SymmetricSparseMatrix<T>& H){
  size_t N = crd.size();
  vector<Matrix3d> diagonal(N);
  for (size_t i = 0; i < N; ++i){
    for (size_t j : cm[i]){
      if (j<=i) continue;
      Vector3d dist = crd[i]-crd[j];
      double d2=normSQ(dist);
      for (int alpha = 0; alpha < 3; ++alpha) {
        double v = SQ(dist[alpha])/d2;
        H.insert(i*3+alpha,j*3+alpha, -v);
        diagonal[i][alpha][alpha]+=v;
        diagonal[j][alpha][alpha]+=v;
        for (int beta = alpha + 1; beta < 3; ++beta) {
          v = dist[alpha]*dist[beta]/d2;
          H.insert(i*3+alpha,j*3+beta, -v);
          H.insert(i*3+beta,j*3+alpha, -v);
          diagonal[i][alpha][beta]+=v;
          diagonal[j][alpha][beta]+=v;
        }
      }
    }
  }

  // insert "diagonal" elements
  for (size_t i = 0; i < N; ++i){
    for (int alpha = 0; alpha < 3; ++alpha) {
      for (int beta = alpha; beta < 3; ++beta) {
        H.insert(i*3+alpha,i*3+beta,diagonal[i][alpha][beta]) ;
      }
    }
  }

}

template<class T>
void integrate_dof(vector<size_t>& intids, ContactMap& newcm, SymmetricSparseMatrix<T>& H){
  // iteratively integrate H'(m,n) = H(m,n) + Sum_k [H(m,k)H(n,k)/H(k,k)]
  for (size_t id : intids){
    for (size_t alpha = 0; alpha < 3; ++alpha) {
      size_t k = id*3+alpha;
      double l_k = 1.0/H(k,k);
      auto& row = H.getRow(k);

      // apply corrections. Note that diagonal elements are recalculated later.
      for (size_t ri = 0; ri < row.size(); ++ri) {
        for (size_t rj = ri+1; rj < row.size(); ++rj) {
          size_t m = row[ri].columnIndex;
          size_t n = row[rj].columnIndex;
          size_t i = m/3;
          size_t j = n/3;
          if (i==j) continue;
          float H_km = row[ri].value;
          float H_kn = row[rj].value;
          float correction = 2*l_k*H_km*H_kn;
          // correct only for i,j in contact
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
      std::vector<SparseMatrixElement<T> >().swap(row);
    }
  }

  // recalculate diagonals for retained beads
  vector<size_t> retainids;
  for (size_t i = 0; i < newcm.size(); ++i)
    if(!is_in(i,intids)) retainids.push_back(i);

  for (size_t i : retainids){
    Matrix3d diagonal;
    for (size_t alpha = 0; alpha < 3; ++alpha) {
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

  // recalculate number of nonzero entries
  size_t nnz = 0;
  for (size_t i = 0; i < H.size(); ++i){
    nnz += H._rows[i].size();
  }
  H._numNonZero = nnz;

}

typedef double real_t;

int main(int argc, char** argv) {
  PdbFile pdb(argv[1]); // already only CA, no check

  vector<Vector3d>& crd = pdb.coords;

  cout << argv[1] << " is a " <<crd.size()<<" atoms file. Computing contact map" <<endl;
  ContactMap cm(crd,7.5);

  SymmetricSparseMatrix<real_t> H0(crd.size()*3);

  cout << "Generating matrix" <<endl;
  generate_matrix(crd,cm,H0);

  H0.saveToFile("interaction_matrix_f.bindat");

  size_t resize_factor = 3;
  vector<size_t> intids, retainids;
  vector<Vector3d> newcrd;
  vector<bool> retain(crd.size(),false);
  for (size_t i = 0; i < crd.size(); ++i){
    if (i%resize_factor == 0) {
      intids.push_back(i);
      newcrd.push_back(crd[i]);
      retain[i]=true;
    }
    else retainids.push_back(i);
  }

  // update cm. In order to save memory, we clear what we don't need
  // and we use the old contact map.

  cout << "Recomputing contact map, resize factor is: " << resize_factor <<endl;
  ContactMap newcm(newcrd,7.5*resize_factor);
  size_t k=0;
  for (size_t i = 0; i < crd.size(); ++i){
    if(retain[i]) {
      cm._cm[i].swap(newcm._cm[k]);
      k++;
    }else {
      std::vector<size_t>().swap(cm._cm[i]);
    }
  }

  cout << "Integrating... " <<endl;
  integrate_dof(intids,cm,H0);
  H0.updateNumNonZero();

  SymmetricSparseMatrix<real_t> H(retainids.size()*3);
  for (size_t i = 0; i < retainids.size(); ++i) {
    for (size_t alpha = 0; alpha < 3; ++alpha) {
      H._rows[i*3+alpha].swap( H0._rows[retainids[i]*3+alpha]);
    }
  }
  H.updateNumNonZero();

  cout << "Done. nnz = " << H.getNumNonZero() <<endl;
  return 0;
}

