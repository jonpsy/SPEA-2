#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief A pseudo code for FineGrainedFitness function
 * Don't compile else thou shalt be buried in errors.
 */
using namespace std;

template <typename MatType>
arma::Col<typename MatType::elem_type>
SPEA2::FineGrainedFitness(const std::vector<arma::Col<ElemType>> &solutionSet)
{
  // Part 1: a) Calculate Strength
  size_t combinedSize = solutionSet.size();
  size_t K = std::ceil(std::sqrt(combinedSize));
  // strength[p] == count of individual p dominates
  arma::uvec strength(combinedSize);
  //dominated[p] == {Indexes in solutionSet who dominate p}
  std::unordered_map<size_t, std::set<size_t> > dominated;

  for(size_t candidateP : len(combinedSize))
  {
    dominated[candidateP] = std::set<size_t>{};

    for(size_t candidateQ : len(combinedSize))
    {
      if(Dominates(solutionSet[candidateP], solutionSet[candidateQ]))
          strength[candidateP]++;
      else if(Dominates(solutionSet[candidateQ], solutionSet[candidateP]))
          dominated[candidateP].insert(candidateQ);
    }
  }

  //Part 1: b) Calculate rawfitness from strength
  arma::vec rawFitness(combinedSize);

  for(size_t candidate : len(combinedSize))
    for(size_t elem : dominated[candidate])
       rawFitness[candidate] += strength[elem];

  arma::vec fineFitness(combinedSize);
  //Part 2: Calculate distance and
  for(size_t idxA = 0; idxA < combinedSize; ++idxA)
  {
    arma::vec distance;
    distance = arma::distance(solutionSet[idxA] - solutionSet.each_col());
    arma::vec sortedDistance = arma::stable_sort(distance);
    fineFitness(idxA) =  rawFitness(idxA) + 1./(sortedDistance(k) + 2.);
  }
}