#include <armadillo>
#include <bits/stdc++.h>

/**
 * @brief A sample code for rawFitness function
 */
using namespace std;

bool Dominates(int a, int b)
{
    return a < b;
}

int main()
{
    size_t populationSize = 4;
    std::vector<int> population(populationSize);
    std::iota(population.begin(), population.end(), populationSize);

    std::map<size_t, size_t> strength; // strength[p] == count of individual p dominates
    std::map<size_t, std::set<size_t> > dominated; //dominated[p] == {Indexes in population who dominate p}

    for(size_t candidateP = 0; candidateP < populationSize; candidateP++)
    {
        strength[candidateP] = 0;
        dominated[candidateP] = std::set<size_t>{};

        for(size_t candidateQ = 0; candidateQ < populationSize; candidateQ++)
        {
            if(Dominates(population[candidateP], population[candidateQ]))
                strength[candidateP]++;
            else if(Dominates(population[candidateQ], population[candidateP]))
                dominated[candidateP].insert(candidateQ);
        }
    }

    std::vector<double> rawFitness(populationSize);

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