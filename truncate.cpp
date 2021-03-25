#include <bits/stdc++.h>
#include <armadillo>

double archiveSize = 10;

// Possibly make it a inner class???
template<typename MatType>
class DistanceSorter
{
 public:
  explicit DistanceSorter(
      const std::vector<std::vector<size_t>> &neighborIndices,
      const std::vector<MatType> &solutionSet)
      : _neighborIndices(neighborIndices),
        _solutionSet(solutionSet) {}

  bool operator()(size_t idxA, size_t idxB)
  {     //..sanity check
    size_t currNeighbor; // The ith neighbor
    double distA, distB;
    // The condition from the paper
    while(currNeighbor < _neighborIndices[0].size() && idxA == idxB)
    {
      distA = arma::distance( //Distance of A to its "i"th neighbor
          _solutionSet[idxA],
          _solutionSet[_neighborIndices[idxA][currNeighbor]]);
      distB = arma::distance( //Distance of B to its "i"th neighbor
          _solutionSet[idxB],
          _solutionSet[_neighborIndices[idxB][currNeighbor]]);
      currNeighbor++;
    }
    // Return idxA if distA < distB
    return distA < distB;
  }

 private:
  std::vector<std::vector<size_t>> &_neighborIndices;
  std::vector<MatType> &_solutionSet;
};


// The size of archive now is numNonDominated
template<typename MatType>
void Truncate(std::vector<MatType>& archive, std::vector<MatType>& solutionSet)
{
  //! Matrix with dynamic size because it re-adjusts
  //! for deleted individuals
  // neighborIndices[i][j] => index of jth closest neighbor of i
  std::vector<std::vector<size_t> > neighborIndices;
  size_t numNonDominated = archive.size();

  for (size_t i : archiveSize)
  { //Cache distance of individual i from others
    arma::vec distances(archiveSize);
    for (size_t j : archiveSize)
      distances[i][j] = arma::norm(archive[i] - archive[j]);
    //Store the indices of the nearest neighbors for individual i
    neighborIndices[i].push_back(arma::stable_sort_index(distances));
  }

  std::vector<size_t> archiveIndices(numNonDominated);
  // This is both the meat of the algorithm and the bottleneck.
  while (archive.size() > archiveSize)
  { // Choose the individual which has the least distance to neighbor.
    size_t idxToRemove = *(std::min_element(
        archiveIndices.begin(), archiveIndices.end(),
        DistanceSorter(neighborIndices, solutionSet)));

    // delete that point from archive and solution set.
    archiveIndices.erase(archiveIndices.begin() + idxToRemove);
    archive.erase(archive.begin() + idxToRemove);
    solutionSet(solutionSet.begin() + idxToRemove);

    // Time to readjust indices
    for (size_t i : neighborIndices.size())
    {
      for (size_t j : neighborIndices.size())
      {
        if (neighborIndices[i][j] == idxToRemove)
        { // Jth neighbor of I is the index to remove
          archive[i].erase(archive.begin() + j);
          solutionSet.erase(solutionSet.begin() + j);
          neighborIndices[i].erase(neighborIndices.begin() + j);
          j--; // To account for the deletion.
        }

        else if (neighborIndices[i][j] > idxToRemove)
        { //Index value decreased by 1 because of deletion.
          neighborIndices[i][j]--;
        }
      }
    }
  }
}