//
// Created by guiyulai on 2021/5/28.
//

#include "problem/mdtlz.h"

#include <cmath>

#include "core/global.h"

namespace emoc {

    mDTLZ1::mDTLZ1(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    mDTLZ1::~mDTLZ1()
    {

    }

    void mDTLZ1::CalObj(Individual *ind)
    {
        int k, size, index;
        double h, g;
        for (int i = 0; i < obj_num_; ++i)
        {
            //h(x_1)
            h = 1.0;
            for (int j = 0; j < obj_num_ - i - 1; ++j)
            {
                h *= ind->dec_[j];
            }
            if (i > 0)
            {
                h *= 1 - ind->dec_[obj_num_ - i - 1];
            }
            h = 0.5 * (1 - h);
            //g(x_2)
            k = i + 1;// [1,m]
            g = 0.0;
            size = (dec_num_ + 1 - k) / obj_num_;
            for (int j = 1; j <= size; ++j)
            {
                index = j * obj_num_ - 1 + k;
                g += pow(ind->dec_[index-1]-0.5,2) - cos(20 * PI * (ind->dec_[index-1]-0.5));
            }
            g = 100 * (size + g);
            ind->obj_[i] = h * (1 + g);
        }
    }

    mDTLZ2::mDTLZ2(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    mDTLZ2::~mDTLZ2()
    {

    }

    void mDTLZ2::CalObj(Individual *ind)
    {
        int k, size, index;
        double h, g;
        for (int i = 0; i < obj_num_; ++i)
        {
            //h(x_1)
            h = 1.0;
            for (int j = 0; j < obj_num_ - i - 1; ++j)
            {
                h *= cos(0.5 * PI * ind->dec_[j]);
            }
            if (i > 0)
            {
                h *= sin(0.5 * PI * ind->dec_[obj_num_ - i - 1]);
            }
            h = 1 - h;
            //g(x_2)
            k = i + 1;// [1,m]
            g = 0.0;
            size = (dec_num_ + 1 - k) / obj_num_;
            for (int j = 1; j <= size; ++j)
            {
                index = j * obj_num_ - 1 + k;
                g += pow(ind->dec_[index-1]-0.5,2);
            }
            ind->obj_[i] = h * (1 + g);
        }
    }

    mDTLZ3::mDTLZ3(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    mDTLZ3::~mDTLZ3()
    {

    }

    void mDTLZ3::CalObj(Individual *ind)
    {
        int k, size, index;
        double h, g;
        for (int i = 0; i < obj_num_; ++i)
        {
            //h(x_1)
            h = 1.0;
            for (int j = 0; j < obj_num_ - i - 1; ++j)
            {
                h *= cos(0.5 * PI * ind->dec_[j]);
            }
            if (i > 0)
            {
                h *= sin(0.5 * PI * ind->dec_[obj_num_ - i - 1]);
            }
            h = 1 - h;
            //g(x_2)
            k = i + 1;// [1,m]
            g = 0.0;
            size = (dec_num_ + 1 - k) / obj_num_;
            for (int j = 1; j <= size; ++j)
            {
                index = j * obj_num_ - 1 + k;
                g += pow(ind->dec_[index-1]-0.5,2) - cos(20 * PI * (ind->dec_[index-1]-0.5));
            }
            g = 100 * (size + g);
            ind->obj_[i] = h * (1 + g);
        }
    }

    mDTLZ4::mDTLZ4(int dec_num, int obj_num) :Problem(dec_num, obj_num)
    {
        for (int i = 0; i < dec_num; ++i)
        {
            lower_bound_[i] = 0.0;
            upper_bound_[i] = 1.0;
        }
    }

    mDTLZ4::~mDTLZ4()
    {

    }

    void mDTLZ4::CalObj(Individual *ind)
    {
        int alpha = 100;
        int k, size, index;
        double h, g;
        for (int i = 0; i < obj_num_; ++i)
        {
            //h(x_1)
            h = 1.0;
            for (int j = 0; j < obj_num_ - i - 1; ++j)
            {
                h *= cos(0.5 * PI * pow(ind->dec_[j], alpha));
            }
            if (i > 0)
            {
                h *= sin(0.5 * PI * pow(ind->dec_[obj_num_ - i - 1], alpha));
            }
            h = 1 - h;
            //g(x_2)
            k = i + 1;// [1,m]
            g = 0.0;
            size = (dec_num_ + 1 - k) / obj_num_;
            for (int j = 1; j <= size; ++j)
            {
                index = j * obj_num_ - 1 + k;
                g += pow(ind->dec_[index-1]-0.5,2);
            }
            ind->obj_[i] = h * (1 + g);
        }
    }
}