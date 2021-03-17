#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief A sample code for Eucledian distance
 * g++ distance.cpp -g -o distance -larmadillo && ./distance
 */

// https://www.youtube.com/watch?v=41e5BIMFD0Q, 31:07 Timestamp, IIT Guwahati NPTEL suggests we should
// also scale (per objective) however pagmo2, and others haven't did that. Further
// the paper doesn't explicitly state that. Can choose to remove it if too complicated.
double scaledDistance(const arma::vec& F1, const arma::vec& F2, const arma::vec& scaleVec)
{
    return std::sqrt(arma::accu(arma::square(F1 - F2) % scaleVec));
}

int main()
{
    //Two vectors in objective space of dimension 2. (Written only for reference, use objectivesMatrix instead)
    arma::vec F1 = {1., 3.};
    arma::vec F2 = {2., 4.};
    arma::mat objectivesMatrix = {{1., 2.}, {3., 4.}}; // stores objective points in column based.
    arma::vec scaleVec(2); //Stores the scale per objective

    scaleVec = 1. / arma::square(arma::max(objectivesMatrix, 1) - arma::min(objectivesMatrix, 1)); //Passed to function

    std::cout << scaledDistance(objectivesMatrix.col(0), objectivesMatrix.col(1), scaleVec) << std::endl;
}