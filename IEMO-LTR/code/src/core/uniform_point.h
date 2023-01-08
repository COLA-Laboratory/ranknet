#pragma once
#include "core/global.h"
#include "core/utility.h"

namespace emoc {

	double** UniformPoint(int num, int *weight_num, int obj_num);
    double** HitAndRun(int num, int obj_num, char *run, char *outputFile);
    double* GoldenPoint(std::string testProblem, int obj_num);
    double** LoadUniformWeights(int num, int obj_num, char* file);

}


