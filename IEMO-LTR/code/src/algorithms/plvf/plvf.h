//
// Created by guiyulai on 2021/6/8.
//

#pragma once
#include "core/individual.h"
#include "algorithms/algorithm.h"
#include "problem/problem.h"

namespace emoc {

    class PLVF : public Algorithm
    {
    public:
        typedef struct
        {
            int index;
            double distance;
        }DistanceInfo;  // store euclidian distance to the index-th weight vector
        typedef struct
        {
            double value;
            int index;
        }SortList;

        PLVF(Problem *problem, int thread_num);
        virtual ~PLVF();

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
        double **rbfnet_c_;
        double *rbfnet_weight_;
        double sigma_;
        double *rbf_output_;
        int mu_;                             // the input size of the rbfnet
        double stepSize_;
    };

}