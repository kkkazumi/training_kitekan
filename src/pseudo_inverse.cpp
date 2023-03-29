#include <stdio.h>
#include <string>
#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

void pseudo_inverse(MatrixXd *O, MatrixXd P){
  MatrixXd Q,R,S;
  Q=P.transpose();
  if( P.rows() > P.cols() ){
    FullPivLU< MatrixXd > ntt(Q*P);
    R=ntt.inverse();
    S=R*Q;
  }
  else{
    FullPivLU< MatrixXd > ntt(P*Q);
    R = ntt.inverse();
    S = Q*R;
  }
  *O=S;
  //使い方
  //pseudo_inverse(&m4,m3);
}

