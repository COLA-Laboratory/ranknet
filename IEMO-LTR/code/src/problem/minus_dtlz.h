#pragma once
#include "problem/problem.h"

namespace emoc {

    class minus_DTLZ1 :public Problem
    {
    public:
        minus_DTLZ1(int dec_num, int obj_num);
        virtual ~minus_DTLZ1();

        void CalObj(Individual *ind);
    };

    class minus_DTLZ2 :public Problem
    {
    public:
        minus_DTLZ2(int dec_num, int obj_num);
        virtual ~minus_DTLZ2();

        void CalObj(Individual *ind);
    };

    class minus_DTLZ3 :public Problem
    {
    public:
        minus_DTLZ3(int dec_num, int obj_num);
        virtual ~minus_DTLZ3();

        void CalObj(Individual *ind);
    };

    class minus_DTLZ4 :public Problem
    {
    public:
        minus_DTLZ4(int dec_num, int obj_num);
        virtual ~minus_DTLZ4();

        void CalObj(Individual *ind);
    };


}