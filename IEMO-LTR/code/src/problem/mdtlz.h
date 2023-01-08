//
// Created by guiyulai on 2021/5/28.
//

#pragma once
#include "problem/problem.h"

namespace emoc {

    class mDTLZ1 :public Problem
    {
    public:
        mDTLZ1(int dec_num, int obj_num);
        virtual ~mDTLZ1();

        void CalObj(Individual *ind);
    };

    class mDTLZ2 :public Problem
    {
    public:
        mDTLZ2(int dec_num, int obj_num);
        virtual ~mDTLZ2();

        void CalObj(Individual *ind);
    };

    class mDTLZ3 :public Problem
    {
    public:
        mDTLZ3(int dec_num, int obj_num);
        virtual ~mDTLZ3();

        void CalObj(Individual *ind);
    };

    class mDTLZ4 :public Problem
    {
    public:
        mDTLZ4(int dec_num, int obj_num);
        virtual ~mDTLZ4();

        void CalObj(Individual *ind);
    };

}
