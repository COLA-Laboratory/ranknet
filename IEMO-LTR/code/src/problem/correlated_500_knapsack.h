//
// Created by gylai on 2021/10/10.
// Contact: Lai Guiyu <guiyulai.chn@gmail.com>
// COLA-Lab@UESTC
//

#ifndef EMOC_CORRELATED_500_KNAPSACK_H
#define EMOC_CORRELATED_500_KNAPSACK_H

#pragma once
#include "core/individual.h"
#include "problem/problem.h"

namespace emoc {
    extern int weight_knapsack[10][500];
    extern int profit_knapsack[10][500];
    extern int capacity_knapsack[10];
    class KNAPSACK:public Problem
    {
    public:
        // note: the smaller the value of alpha, the more correlated the objectives (with f1 and f2) are.
        KNAPSACK(int dec_num, int obj_num, double alpha = 0.1);
        virtual ~KNAPSACK();

        void CalObj(Individual *ind);

    public:
        double alpha;// correlated parameter
    };

}

#endif //EMOC_CORRELATED_500_KNAPSACK_H
