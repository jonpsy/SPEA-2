#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief Polynomial mutation
 * Code length has been decreased significantly.
 * Inspired from DEAP/mutation.py https://github.com/DEAP/deap/blob/d328fe6b68e7528b2d2d990bb2ab1ad1786e6f58/deap/tools/mutation.py#L51
 * Naming convention followed.
 * g++ mutation.cpp -g -o mutation -larmadillo
 *
 */

//Variables to be defined for polynomial mutation
double perturbationIndex = 30; //eta_m
double mutationProb = 1; // For testing
size_t numVariables = 5;
arma::vec lowerBound {-10, -10, -10, -10, -10};
arma::vec upperBound {10, 10, 10, 10, 10};

void enforceBounds(arma::vec& child)
{
  child = arma::min(arma::max(child, lowerBound), upperBound);
}

void PolynomialMutation(const arma::vec& parent, arma::vec& child)
{
  for(size_t geneIdx=0; geneIdx < numVariables; ++geneIdx)
  {
    if(arma::randu() > mutationProb) //Should these gene be mutated?
      continue;

    double geneRange = upperBound[geneIdx] - lowerBound[geneIdx];
    //Normalised distance from bounds
    double lowerDelta = (parent[geneIdx] - lowerBound[geneIdx]) / geneRange;
    double upperDelta = (upperBound[geneIdx] - parent[geneIdx]) / geneRange;
    double mutationPower = 1. / (perturbationIndex + 1.0);
    double value; //Intermediatery
    double perturbationFactor;
    double rand = arma::randu();
    if(rand < 0.5)
    {
      value = 2. * rand + (1. - 2. * rand) *
          std::pow(upperDelta, perturbationIndex + 1.0);
      perturbationFactor = std::pow(value, mutationPower) - 1.;
    }
    else
    {
      value = 2. * (1. - rand) + 2.*(rand - 0.5) *
          std::pow(lowerDelta, perturbationIndex + 1.0);
      perturbationFactor = 1. - std::pow(value, mutationPower);
    }

    child[geneIdx] += perturbationFactor * geneRange;
  }
  enforceBounds(child);
}

//Vectorised WIP!
void PolynomialMutationVector(const arma::vec& parent, arma::vec& child)
{
  arma::vec geneRange = upperBound- lowerBound;
  arma::vec lowerDelta = (parent - lowerBound) / geneRange;
  arma::vec upperDelta = (upperBound - parent) / geneRange;
  double mutationPower = 1. / (perturbationIndex + 1.0);
  arma::vec value(numVariables, arma::fill::zeros);
  arma::vec perturbationFactor;
  arma::vec rand = arma::randu(numVariables, 1);
  value = 2. * rand + (1. - 2. * rand) *
          arma::pow(upperDelta, 1 / mutationPower) % (rand < 0.5);
  value += 2. * (1. - rand) + 2.*(rand - 0.5) *
          arma::pow(lowerDelta, 1 / mutationPower) % (rand >= 0.5);

  // perturbationFactor arma::pow(value, mutationPower) - 1.
  child += perturbationFactor* geneRange;
  // chidl %= arma::randu
  enforceBounds(child);
}

int main()
{
  arma::vec parentA { 1., 2., 3., 4., 5.};
  arma::vec childA = parentA;
  PolynomialMutation(parentA, childA);

  std::cout << childA << std::endl;
}