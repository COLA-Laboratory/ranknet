//
// Created by guiyulai on 2021/7/19.
//
#pragma once
#include <vector>
#include <numeric>
#include <queue>
#include <algorithm>

#include "core/individual.h"
#include "algorithms/algorithm.h"
#include "problem/problem.h"
#include "core/global.h"
#include "core/nd_sort.h"
#include "core/tournament_selection.h"
#include "operator/mutation.h"
#include "operator/sbx.h"
#include "random/random.h"
#include "core/utility.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>

namespace emoc {

    class NSGA2_LTR : public Algorithm
    {
    public:
        typedef struct
        {
            int index;
            double distance;
        }DistanceInfo;  // store crowding distance of index-th individual
        typedef struct
        {
            int index;
            double value;
        }SortList;

        NSGA2_LTR(Problem *problem, int thread_num);
        virtual ~NSGA2_LTR();

        void Run();

    private:
        void Initialization();
        void Crossover(Individual **parent_pop, Individual **offspring_pop);

        // set the crowding distance of given individual index
        void SetDistanceInfo(std::vector<DistanceInfo> &distanceinfo_vec, int target_index, double distance);
        void SetPreferenceInfo(std::vector<SortList> &preferenceInfo_vec, int target_index, double value);

        // use crowding distance to sort the individuals with rank rank_index, the result is stored in pop_sort
        // return the number of indviduals of rank rank_index
        int CrowdingDistance(Individual **mixed_pop, int pop_num, int *pop_sort, int rank_index);
        /**
        * select arbitrarily from the last PF
        * @param mixed_pop
        * @param pop_num
        * @param pop_sort
        * @param rank_index
        * @return the number of indviduals of rank rank_index
        */
        int ArbitrarilySelect(Individual **mixed_pop, int pop_num, int *pop_sort, int rank_index);
        /**
         * use LP to sort the individuals with rank rank_index, the result is stored in pop_sort
         * @param mixed_pop
         * @param pop_num
         * @param pop_sort
         * @param rank_index
         * @return the number of indviduals of rank rank_index
         */
        int RankViaPreference(Individual **mixed_pop, int pop_num, int *pop_sort, int rank_index, int nneed);
        // do nsga2's environment selection on mixed_pop, the result is stored in parent_pop
        bool EnvironmentalSelection(Individual **parent_pop, Individual **mixed_pop, bool usingPreferenceModel);

        // solve lp
    private:
        double *weight_;                   // weight for the DM
        int weight_num_;                   // the number of weight vector
        std::queue<double *> winners_;
        std::queue<double *> losers_;
    };

}
