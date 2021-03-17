#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief A sample code for rawFitness function
 * g++ raw_fit.cpp -g -o raw_fit -larmadillo && ./raw_fit
 */
using namespace std;

bool Dominates(int a, int b)
{
    return a < b;
}

int main()
{
    size_t populationSize = 4;
    int startPoint = 3;
    std::vector<int> population(populationSize);
    std::iota(population.begin(), population.end(), startPoint);

    arma::uvec strength(populationSize, arma::fill::zeros); // strength[p] == count of individual p dominates
    std::map<size_t, std::set<size_t> > dominated; //dominated[p] == {Indexes in population who dominate p}

    for(size_t candidateP = 0; candidateP < populationSize; candidateP++)
    {
        dominated[candidateP] = std::set<size_t>{};

        for(size_t candidateQ = 0; candidateQ < populationSize; candidateQ++)
        {
            if(Dominates(population[candidateP], population[candidateQ]))
                strength[candidateP]++;
            else if(Dominates(population[candidateQ], population[candidateP]))
                dominated[candidateP].insert(candidateQ);
        }
    }

    arma::vec rawFitness(populationSize, arma::fill::zeros);

    for(size_t candidate =0; candidate < populationSize; candidate++)
    {
        for(const size_t& elem : dominated[candidate])
        {
            rawFitness[candidate] += strength[elem];
        }
    }

    for(const auto& elem : rawFitness)
    {
        cout << elem << endl;
    }
}