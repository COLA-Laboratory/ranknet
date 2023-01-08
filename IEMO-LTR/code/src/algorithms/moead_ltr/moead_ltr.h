//
// Created by guiyulai on 2021/6/9.
//

#pragma once

#define INCREASING_STEP_SIZE 0 // whether to increase the step size gradually
#include <vector>
#include <numeric>
#include <queue>

#include "core/individual.h"
#include "algorithms/algorithm.h"
#include "problem/problem.h"
#include "core/global.h"

namespace emoc {

    class MOEAD_LTR : public Algorithm
    {
    public:
        typedef struct
        {
            int index;
            double distance;
        }DistanceInfo;  // store euclidian distance to the index-th weight vector

        MOEAD_LTR(Problem *problem, int thread_num);
        virtual ~MOEAD_LTR();

        void Run();

    private:
        void Initialization();
        void InitializeNeighbours();
        void UpdateNeighbours();
        void Crossover(Individual **parent_pop, int current_index, Individual *offspring);

        // use offspring to update the neighbour of current_index-th individual with specified aggregation function
        void UpdateSubproblem(Individual *offspring, int current_index, int aggregation_type);

        void Consultation_PreferenceElicitation();
        void UsingGoldenFunc();
        void RecordCurrentPop(PyObject *pop);
        void UpdateTrainingSet();
        void LoadTrainingSet(PyObject *winners, PyObject *losers);
        PyObject *TrainRankNet_ReturnScore(PyObject *pFunction,PyObject *winners, PyObject *losers, PyObject *currPop);
        void UpdateScoreByRankNet(PyObject * res, double *score);
        void BiasingWeightSet();


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
        double step_size_;                  // the speed of moving the weight vectors
        double retention_rate_;             // proportion of promising weigth vectors
        std::queue<double *> winners_;
        std::queue<double *> losers_;
        double *score_RankNet_;
    };

}
