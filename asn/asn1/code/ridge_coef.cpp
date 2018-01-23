#include <Rcpp.h>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
using namespace Rcpp;
using namespace Eigen;

template <class out,class inp>
inline out convertMatrix(const inp& matrixinput) {
  int rows = matrixinput.rows();
  int cols = matrixinput.cols();
  out matrixOutput(rows, cols);
  for(int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      matrixOutput(i, j) = matrixinput(i, j);
    }
  }
  return matrixOutput;
}
template <class out, class inp>
inline out convertVector(const inp& vectorinput) {
  int len = vectorinput.size();
  out vectoroutput(len);
  for (int i = 0; i < len; ++i) {
    vectoroutput(i) = vectorinput(i);
  }
  return vectoroutput;
}

// [[Rcpp::export]]
NumericMatrix ridge_coef_cpp(NumericMatrix X, NumericVector y, double l) {
  MatrixXd A = convertMatrix<Eigen::MatrixXd, Rcpp::NumericMatrix>(X);
  VectorXd b = convertVector<Eigen::VectorXd, Rcpp::NumericVector>(y);
  MatrixXd I = VectorXd::Ones(A.cols()).asDiagonal();
  
  MatrixXd At = A.transpose();
  LLT<MatrixXd> llt;
  llt.compute(At * A + l * I);
  
  return convertMatrix<Rcpp::NumericMatrix, Eigen::MatrixXd>(llt.solve(At * b));
}
