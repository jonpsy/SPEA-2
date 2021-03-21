#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief Simulated Binary CrossOver.
 * Paper: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.7291&rep=rep1&type=pdf
 * Appendix A : SBX
 * Inspired from NSGA-II Code by Deb in C
 * Gist: https://gist.github.com/Tiagoperes/1779d5f1c89bae0cfdb87b1960bba36d by @Tiagoperes
 *
 * g++ sbx.cpp -g -o sbx -larmadillo && ./sbx
 */

// Variables global here are defined in class ctor or loop
double crossoverProb = 1; // delta
double spreadIndex = 2; // eta_c
double epsilon = 1e-14;
double numVariables = 5;
arma::vec lowerBound {-10, -10, -10, -10, -10};
arma::vec upperBound {10, 10, 10, 10, 10};

void enforceBounds(double& child, double lower, double upper)
{
  child = std::min(std::max(child, lower), upper);
}

void SimulatedBinaryCrossover(const arma::vec& parentA,const arma::vec& parentB, arma::vec& childA, arma::vec& childB)
{

  if (arma::randu() > crossoverProb) //Should I crossover?
    return;

  for (size_t geneIdx=0; geneIdx < numVariables; ++geneIdx)
  {
    if (arma::randu() > 0.5) // Each gene has 50% of changing its value
      continue;
    if(std::fabs(parentA[geneIdx] - parentB[geneIdx]) < epsilon) // If parents are equal, skip
      continue;

    double recessiveGene = std::min(parentA[geneIdx], parentB[geneIdx]);
    double dominantGene = std::max(parentA[geneIdx], parentB[geneIdx]);

    double geneRange = dominantGene - recessiveGene;
    double geneAverage = 0.5* (dominantGene + recessiveGene);
    //! Calculate gene for childA
    //! a) Variables for childA (recessive - lowerbound)
    double beta = 1. + 2. * (recessiveGene - lowerBound[geneIdx])/ geneRange;
    double alpha = 2. - std::pow(beta, -(spreadIndex + 1.));
    double spreadFactor;
    double rand = arma::randu();
    //! b) spreadFactor for bounded SBX
    if (rand <= 1./alpha)
    {
      spreadFactor = std::pow(alpha * rand, 1./(spreadIndex + 1.));
    }
    else
    {
      spreadFactor = std::pow(1./(2. - rand * alpha), 1./(spreadIndex + 1.));
    }

    //! c) Fill child A
    childA[geneIdx] = geneAverage - 0.5 * spreadFactor * geneRange;

    //! Calculate for childB
    //! a) Variables for childB (upperbound - dominant)
    beta = 1. + 2. * (upperBound[geneIdx] - dominantGene)/ geneRange;
    alpha = 2. + std::pow(beta, -(spreadIndex + 1.));
    //! b) spreadFactor for bounded SBX
    if (rand <= 1./alpha)
    {
      spreadFactor = std::pow(alpha * rand, 1./(spreadIndex + 1.));
    }
    else
    {
      spreadFactor = std::pow(1./(2. - rand * alpha), 1./(spreadIndex + 1.));
    }

    //! c)Fill childB
    childB[geneIdx] = geneAverage  + 0.5 * spreadFactor * geneRange;

    //Bound check
    enforceBounds(childA[geneIdx], lowerBound[geneIdx], upperBound[geneIdx]);
    enforceBounds(childB[geneIdx], lowerBound[geneIdx], upperBound[geneIdx]);
  }
}

//! WIP
void SBX_Vectorised(arma::vec parentA, arma::vec parentB, arma::vec& childA, arma::vec& childB)
{
  arma::uvec idxEqual = (parentA == parentB);
  // arma::uvec idxProbability = (arma::randu(2, 1) < 0.5);
  arma::vec recessiveGene = arma::min(parentA, parentB);
  arma::vec dominantGene = arma::max(parentA, parentB);
  arma::vec geneAverage = 0.5 * (parentA + parentB);
  arma::vec geneRange = (dominantGene - recessiveGene);   //Some value becomes nan (dominanteGene[i] == recessiveGene[i])
  arma::vec beta = (1. + 2. * (recessiveGene - lowerBound)/ (dominantGene - recessiveGene));
  arma::vec alpha = 2. - arma::pow(beta, - (spreadIndex + 1.));
  arma::vec rand = arma::randu(2, 1);
  arma::vec spreadFactor = arma::pow(alpha * rand, 1./(spreadIndex + 1.));
  childA = geneAverage + 0.5 * spreadFactor * geneRange;
  childA = childA % (1 - idxEqual) + parentA % (1 - idxEqual); //nan replaced with parent
  // arma::vec spreadFactor = arma::pow(alpha * rand)
  // arma::uvec idxSpreadFactor = arma::randu(2, 1) < 1./alpha;

}

int main()
{
  arma::vec parentA {1., 2. ,3., 4., 5., 6.};
  arma::vec parentB {1., 4., 7., 8., 9., 10.};
  arma::vec childA = parentA;
  arma::vec childB = parentB;

  SimulatedBinaryCrossover(parentA, parentB, childA, childB);

  std::cout << "CHILD A: " << childA << std::endl;
  std::cout << "CHILD B: " << childB << std::endl;
  // std::cout <<  arma::randu(2, 1) << std::endl;
}