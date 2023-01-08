//
// Created by guiyulai on 2021/5/28.
//

#include "problem/minus_dtlz.h"

#include <cmath>

#include "core/global.h"

namespace emoc {

    minus_DTLZ1::minus_DTLZ1(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    minus_DTLZ1::~minus_DTLZ1()
    {

    }

    void minus_DTLZ1::CalObj(Individual *ind)
    {
        double gx = 0.0;
        int k = dec_num_ - obj_num_ + 1;

        for (int i = dec_num_ - k; i < dec_num_; i++)
            gx += pow((ind->dec_[i] - 0.5), 2.0) - cos(20.0 * PI * (ind->dec_[i] - 0.5));
        gx = 100.0 * (k + gx);

        for (int i = 0; i < obj_num_; i++)
            ind->obj_[i] = 0.5 * (1.0 + gx);

        for (int i = 0; i < obj_num_; i++)
        {
            for (int j = 0; j < obj_num_ - (i + 1); j++)
                ind->obj_[i] *= ind->dec_[j];
            if (i != 0)
            {
                int aux = obj_num_ - (i + 1);
                ind->obj_[i] *= 1 - ind->dec_[aux];
            }
            ind->obj_[i] = -1.0 * ind->obj_[i];
        }
    }

    minus_DTLZ2::minus_DTLZ2(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    minus_DTLZ2::~minus_DTLZ2()
    {

    }

    void minus_DTLZ2::CalObj(Individual *ind)
    {
        double gx = 0.0;
        int k = dec_num_ - obj_num_ + 1;

        for (int i = dec_num_ - k; i < dec_num_; i++)
            gx += pow((ind->dec_[i] - 0.5), 2.0);

        for (int i = 0; i < obj_num_; i++)
            ind->obj_[i] = -1.0 * (1.0 + gx);

        for (int i = 0; i < obj_num_; i++)
        {
            for (int j = 0; j < obj_num_ - (i + 1); j++)
                ind->obj_[i] *= cos(PI * 0.5 * ind->dec_[j]);
            if (i != 0)
            {
                int aux = obj_num_ - (i + 1);
                ind->obj_[i] *= sin(PI * 0.5 * ind->dec_[aux]);
            }
        }
    }

    minus_DTLZ3::minus_DTLZ3(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    minus_DTLZ3::~minus_DTLZ3()
    {

    }

    void minus_DTLZ3::CalObj(Individual *ind)
    {
        double gx = 0.0;
        int k =dec_num_ -obj_num_ + 1;

        for (int i =dec_num_ - k; i <dec_num_; i++)
            gx += pow((ind->dec_[i] - 0.5), 2.0) - cos(20.0 * PI * (ind->dec_[i] - 0.5));
        gx = 100.0 * (k + gx);

        for (int i = 0; i <obj_num_; i++)
            ind->obj_[i] = -1.0 * (1.0 + gx);

        for (int i = 0; i <obj_num_; i++)
        {
            for (int j = 0; j <obj_num_ - (i + 1); j++)
                ind->obj_[i] *= cos(PI * 0.5 * ind->dec_[j]);
            if (i != 0)
            {
                int aux =obj_num_ - (i + 1);
                ind->obj_[i] *= sin(PI * 0.5 * ind->dec_[aux]);
            }
        }
    }

    minus_DTLZ4::minus_DTLZ4(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    minus_DTLZ4::~minus_DTLZ4()
    {

    }

    void minus_DTLZ4::CalObj(Individual *ind)
    {
        double gx = 0.0, alpha = 100.0;
        int k =dec_num_ -obj_num_ + 1;

        for (int i =dec_num_ - k; i <dec_num_; i++)
            gx += pow((ind->dec_[i] - 0.5), 2.0);

        for (int i = 0; i <obj_num_; i++)
            ind->obj_[i] = -1.0 * (1.0 + gx);

        for (int i = 0; i <obj_num_; i++)
        {
            for (int j = 0; j <obj_num_ - (i + 1); j++)
                ind->obj_[i] *= cos(PI * 0.5 * pow(ind->dec_[j], alpha));
            if (i != 0)
            {
                int aux =obj_num_ - (i + 1);
                ind->obj_[i] *= sin(PI * 0.5 * pow(ind->dec_[aux], alpha));
            }
        }
    }

}
