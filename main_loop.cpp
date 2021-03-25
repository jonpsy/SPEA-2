  #include <bits/stdc++.h>
  #include <armadillo>

  template<typename MatType,
          typename... ArbitraryFunctionType,
          typename... CallbackTypes>
  typename MatType::elem_type SPEA2::Optimize(
      std::tuple<ArbitraryFunctionType...>& objectives,
      MatType& iterate,
      CallbackTypes&&... callbacks)
  {
    // .. Sanity checks & typedef
    typedef std::reference_wrapper<const MatType> ConstRefMatType;
    size_t combinedSize = populationSize + archiveSize;
    // Archive size is same at the end of each iteration.
    std::vector<MatType> population(populationSize); //P_{0}
    std::vector<MatType> archive(archiveSize);// A_{0} = empty
    // Modifying operations are never done on combinedPopulation
    // Better to store by const ref. population and archive :)
    // C_{0}=P_{0} U A_{0}
    std::vector<ConstRefMatType> combinedPopulation(combinedSize);
    std::vector<arma::Col<ElemType>> solutionSet(combinedSize);
    // Initialize uniformly around the starting point
    population = arma::randu(cols, rows) - 0.5 + iterate;
    // We dont store archive because its empty
    StorePopulation(combinedPopulation, population);
    // Store the corresponding objective values to solution set.
    EvaluateObjectives(combinedPopulation, objectives,
        solutionSet);

    // Iteratively upgrade the archive population and return it.
    while (gen = 0; gen <= maxGen && terminate != true; ++gen)
    {       // Loop runs for maxGen + 1 times
      terminate |= Callback::StepTaken(...);

      //! [1] Fitness Assignment.
      // FineGrainedFitness => A two step fitness calculation method.
      arma::Col<ElemType> objectiveFitness =
          FineGrainedFitness(solutionSet);
      // Lower fitness is better.
      arma::uvec sortedIndices = arma::stable_sort_index(objectiveFitness);
      // num solutions for which objectiveFitness is < 1.
      size_t numNonDominated = std::count(objectiveFitness,
          [&](Elemtype a){ a < 1; })

      //! [2] Environment Selection.
      // Copy non dominated solutions to archive.
      if (numNonDominated > archiveSize) // Diversity preserve.
        archive.clear();
        std::copy(archive, combinedPopulation);
        // Truncate so that numNonDominated == archive
        Truncate(archive);
      else //Fill in order of descending fitness.
        archive = combinedPopulation[sortedIndices(
            arma::span(0, archiveSize))];

      EvaluateObjectives(combinedPopulation, objectives,
          solutionSet);

      // No further modification, just return archive A_{t+1}
      if (generation == maxGeneration)
        break;

      // Genetic operation to get P_{t+1}
      ModifyPopulation(population, archive, objectiveFitness);
      StorePopulation(combinedPopulation, population);
      StoreArchive(combinedPopulation, archive);
    }

    bestFront = archive;
    CallBack::EndOptimization(...);
    // Store the minimum sum of objectives.
    performance = arma::accu(min_accumulation(solutionSet));

    return performance;
  }