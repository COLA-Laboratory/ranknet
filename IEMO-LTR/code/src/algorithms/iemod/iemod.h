//
// Created by guiyulai on 2021/5/17.
//

#pragma once
#include "core/individual.h"
#include "algorithms/algorithm.h"
#include "problem/problem.h"

namespace emoc {

    class IEMOD : public Algorithm
    {
    public:
        typedef struct
        {
            int index;
            double distance;
        }DistanceInfo;  // store euclidian distance to the index-th weight vector

        IEMOD(Problem *problem, int thread_num);
        virtual ~IEMOD();

        void Run();

    private:
        void Initialization();
        void InitializeNeighbours();// initialize neighborhood, including allocate the memory
        void UpdateNeighbours();
        void Crossover(Individual **parent_pop, int current_index, Individual *offspring);

        // use offspring to update the neighbour of current_index-th individual with specified aggregation function
        void UpdateSubproblem(Individual *offspring, int current_index, int aggregation_type);


    private:
        double **lambda_;                  // weight vector
        double *weight_;                   // weight for the DM
        double *goldenPoint_;              // the golden point
        int weight_num_;                   // the number of weight vector
        int **neighbour_;	               // neighbours of each individual
        int neighbour_num_;                // the number of neighbours
        double *ideal_point_;
        int aggregation_type_;
        double pbi_theta_;
        double alpha_;
        double alpha_DM_;
    };

}
