#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief A sample code for Eucledian distance
 * g++ distance.cpp -g -o distance -larmadillo && ./distance
 */

double printDistance(const arma::vec& F1, const arma::vec& F2, const arma::vec& scaleVec)
{
    return std::sqrt(arma::accu(arma::square(F1 - F2) % scaleVec));
}

int main()
{
    //Two vectors in objective space of dimension 2. (Written only for reference, use fMatrix instead)
    arma::vec F1 = {1., 3.};
    arma::vec F2 = {2., 4.};
    arma::mat fMatrix = {{1., 2.}, {3., 4.}}; // stores objective points in column based.
    arma::vec scaleVec(2); //Stores the scale per objective

    scaleVec = 1. / arma::square(arma::max(fMatrix, 1) - arma::min(fMatrix, 1)); //Passed to function

    std::cout << printDistance(fMatrix.col(0), fMatrix.col(1), scaleVec) << std::endl;
}