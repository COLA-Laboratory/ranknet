//
// Created by gylai on 2021/7/27.
//

#pragma once
#include <vector>
#include <numeric>
#include <queue>
#include "core/individual.h"
#include "algorithms/algorithm.h"
#include "problem/problem.h"
#include "core/global.h"


namespace emoc {

    class R2_IBEA_LTR : public Algorithm
    {
    public:
        R2_IBEA_LTR(Problem *problem, int thread_num);
        virtual ~R2_IBEA_LTR();

        void Run();

    private:
        void Initialization();
        void Crossover(Individual **parent_pop, Individual **offspring_pop);

        double CalEpsIndicator(Individual *ind1, Individual *ind2);
        void CalFitness(Individual **pop, int pop_num, double *fitness);
        void EnvironmentalSelection(Individual **parent_pop, Individual **mixed_pop);

        double CalR2Indicator(Individual *x);
        double CalR2Indicator(Individual *x, Individual *y);
        double CalBiR2Indicator(Individual *x, Individual *y);
        Individual* TournamentByBiR2Indicator(Individual *ind1, Individual *ind2);
        void UpdateTrainingSet();
        void RecordCurrentPop(PyObject *pop);
        void LoadTrainingSet(PyObject *winners, PyObject *losers);
        PyObject *TrainRankNet_ReturnScore(PyObject *pFunction,PyObject *winners, PyObject *losers, PyObject *currPop);
        void UpdateScoreByRankNet(PyObject * res, double *score);
        void BiasingWeightSet();
        void Consultation_PreferenceElicitation();

    private:
        double kappa_;
        double *reference_point_;
        double **lambda_;
        double *weight_;
        double step_size_;
        double *ideal_point_;
        double retention_rate_;
        std::queue<double *> winners_;
        std::queue<double *> losers_;
        double *score_RankNet_;
    };

}
