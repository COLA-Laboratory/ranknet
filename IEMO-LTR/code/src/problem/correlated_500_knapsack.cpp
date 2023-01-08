//
// Created by gylai on 2021/10/10.
// Contact: Lai Guiyu <guiyulai.chn@gmail.com>
// COLA-Lab@UESTC
//

#include "correlated_500_knapsack.h"
#include <cmath>
#include "core/global.h"

namespace emoc {
    int weight_knapsack[10][500];// TODO: 10 --> dynamic
    int profit_knapsack[10][500];
    int capacity_knapsack[10];
    KNAPSACK::KNAPSACK(int dec_num, int obj_num, double alpha)
    :Problem(dec_num, obj_num),
    alpha(alpha)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    KNAPSACK::~KNAPSACK()
    {

    }

    void KNAPSACK::CalObj(Individual *ind)
    {
        // for bi-dec_[i], the value is 0 or 1
        assert(obj_num_<=10);
        double random_f[10]={0};
        for (int i = 0; i < obj_num_; ++i)
        {
            for (int j = 0; j < 500; ++j)
            {
                random_f[i] += ind->dec_[j] * profit_knapsack[i][j];
            }
        }
        for (int i = 0; i < obj_num_; ++i)
        {
            if (i < 2)
            {
                ind->obj_[i] = random_f[i];
            }else if (i % 2 == 0)
            {
                ind->obj_[i] = alpha * random_f[i] + (1-alpha) * random_f[0];
            } else{
                ind->obj_[i] = alpha * random_f[i] + (1-alpha) * random_f[1];
            }
        }
        // maximize problem --> minimize problem
        for (int i = 0; i < obj_num_; ++i)
        {
            ind->obj_[i] = -1.0 * ind->obj_[i];
        }
    }
}