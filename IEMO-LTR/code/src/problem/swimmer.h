//
// Created by gylai on 22-9-28.
// Contact: Lai Guiyu <guiyulai.chn@gmail.com>
// COLA-Lab@UESTC
//

#pragma once
#include "core/individual.h"
#include "problem/problem.h"

#define PY_SSIZE_T_CLEAN
#include <Python.h>

namespace emoc {

    class SWIMMER:public Problem
    {
    public:
        SWIMMER(int dec_num, int obj_num);
        virtual ~SWIMMER();

        void CalObj(Individual *ind);
    };

}