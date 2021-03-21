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
double crossoverProb = 0.001; //Keeping it very less for testing
double epsilon = 1e-14;
double distributionIndex = 2; // eta_c
double numVariables = 2;
arma::vec lowerBound {-10, -10};
arma::vec upperBound {10, 10};

void enforceBounds(double& child, double lower, double upper)
{
  child = std::min(std::max(child, lower), upper);
}

void SimulatedBinaryCrossover(const arma::vec& parentA,const arma::vec& parentB, arma::vec& childA, arma::vec& childB)
{
  double recessiveGene;
  double dominantGene;

  if (arma::randu() < crossoverProb) //Should I crossover? For testing, definetly!
  {
    for (size_t geneIdx=0; geneIdx < numVariables; ++geneIdx)
    {
      if (arma::randu() <= 0.5) // 50% chance of gene changing value
      {
        if(std::fabs(parentA[geneIdx] - parentB[geneIdx]) > epsilon)
        {

          recessiveGene = std::min(parentA[geneIdx], parentB[geneIdx]);
          dominantGene = std::max(parentA[geneIdx], parentB[geneIdx]);

          double geneRange = dominantGene - recessiveGene;
          double geneAverage = 0.5* (dominantGene + recessiveGene);
          // Calculate child A
          double beta = 1. + 2. * (recessiveGene - lowerBound[geneIdx])/ geneRange;
          double alpha = 2. - std::pow(beta, -(distributionIndex + 1.));
          double spreadFactor;
          double rand = arma::randu();
          if (rand <= 1./alpha)
          {
            spreadFactor = std::pow(alpha * rand, 1./(distributionIndex + 1.));
          }
          else
          {
            spreadFactor = std::pow(1./(2. - rand * alpha), 1./(distributionIndex + 1.));
          }

          // Fill child A
          childA[geneIdx] = geneAverage - 0.5 * spreadFactor * geneRange;

          //Calculate for childB
          beta = 1. + 2. * (upperBound[geneIdx] - dominantGene)/ geneRange;
          alpha = 2. + std::pow(beta, -(distributionIndex + 1.));

          if (rand <= 1./alpha)
          {
            spreadFactor = std::pow(alpha * rand, 1./(distributionIndex + 1.));
          }
          else
          {
            spreadFactor = std::pow(1./(2. - rand * alpha), 1./(distributionIndex + 1.));
          }

          //Fill childB
          childB[geneIdx] = geneAverage  + 0.5 * spreadFactor * geneRange;

          //Bound check
          enforceBounds(childA[geneIdx], lowerBound[geneIdx], upperBound[geneIdx]);
          enforceBounds(childB[geneIdx], lowerBound[geneIdx], upperBound[geneIdx]);
        }
      }
    }
  }
}

//! WIP
void SBX_Vectorised(arma::vec parentA, arma::vec parentB, arma::vec& childA, arma::vec& childB)
{
  arma::uvec idxEqual = (parentA == parentB);
  arma::uvec idxProbability = (arma::randu(2, 1) < 0.5);
 // Get genes from parentA where less + genes from when B less, when equal remains zero
  arma::vec recessiveGene = (parentA < parentB) % parentA + (parentB < parentA) % parentB;
  arma::vec dominantGene = (parentA > parentB) % parentA + (parentB > parentA) % parentB;
  arma::vec geneRange = (dominantGene - recessiveGene);   //Some value becomes nan (dominanteGene[i] == recessiveGene[i])
  arma::vec beta = (1. + 2. * (recessiveGene - lowerBound)/ (dominantGene - recessiveGene));
  arma::vec alpha = 2. - arma::pow(beta, - (distributionIndex + 1.));
  arma::uvec idxSpreadFactor = arma::randu(2, 1) < 1./alpha;
}

int main()
{
  arma::vec parentA {1., 2.};
  arma::vec parentB {1., 4.};
  arma::vec childA = parentA;
  arma::vec childB = parentB;

  SimulatedBinaryCrossover(parentA, parentB, childA, childB);

  std::cout << "CHILD A: " << childA << std::endl;
  std::cout << "CHILD B: " << childB << std::endl;
  std::cout <<  childA << std::endl;
}