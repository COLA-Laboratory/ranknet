//
// Created by gylai on 2021/7/27.
//

#include "algorithms/r2_ibea_ltr/r2_ibea_ltr.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>

#include "core/global.h"
#include "core/tournament_selection.h"
#include "operator/sbx.h"
#include "operator/mutation.h"
#include "random/random.h"

#include "core/utility.h"
#include "core/uniform_point.h"


namespace emoc {


    R2_IBEA_LTR::R2_IBEA_LTR(Problem *problem, int thread_num):
    Algorithm(problem, thread_num),
    kappa_(0.005),
    ideal_point_(new double[g_GlobalSettings->obj_num_]),
    reference_point_(new double[g_GlobalSettings->obj_num_]),
    lambda_(nullptr),
    weight_(nullptr),
    step_size_(0.2),
    retention_rate_(0.2),
    score_RankNet_(new double[g_GlobalSettings->population_num_])
    {

    }

    R2_IBEA_LTR::~R2_IBEA_LTR()
    {
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            delete[] lambda_[i];
            lambda_[i] = nullptr;
        }
        delete[] lambda_;
        delete[] weight_;
        delete[] ideal_point_;
        delete[] reference_point_;
        delete[] score_RankNet_;
        lambda_ = nullptr;
        weight_ = nullptr;
        ideal_point_ = nullptr;
        reference_point_ = nullptr;
        score_RankNet_ = nullptr;

        while(winners_.size())
        {
            delete[] winners_.front();
            winners_.pop();
            delete[] losers_.front();
            losers_.pop();
        }
    }

    void R2_IBEA_LTR::Run()
    {
        // TODO: multi-thread for Python
        class PyThreadStateLock PyThreadLock;

        Initialization();
        if (g_GlobalSettings->output_interval_)
        {
            TrackPopulation(g_GlobalSettings->iteration_num_);
        }
        while (!g_GlobalSettings->IsTermination())
        {
            // begin each iteration
            g_GlobalSettings->iteration_num_++;

            // consultation, preference elicitation
            if (g_GlobalSettings->tau_&&(g_GlobalSettings->iteration_num_ % g_GlobalSettings->tau_ == 0))
            {
                Consultation_PreferenceElicitation();
            }
            /**Py_BEGIN_ALLOW_THREADS, Py_END_ALLOW_THREADS
             * 为了在较长时间的C函数调用前，
             * 临时释放全局锁，
             * 完成后重新获取全局锁，
             * 以避免阻塞其他python的线程继续运行
             */
            // Py_BEGIN_ALLOW_THREADS;

            // generate offspring population
            Crossover(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->offspring_population_.data());
            MutationPop(g_GlobalSettings->offspring_population_.data(), 2 * g_GlobalSettings->population_num_ / 2, g_GlobalSettings);
            EvaluatePop(g_GlobalSettings->offspring_population_.data(), 2 * g_GlobalSettings->population_num_ / 2);
            MergePopulation(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_, g_GlobalSettings->offspring_population_.data(),
                            2 * g_GlobalSettings->population_num_ / 2, g_GlobalSettings->mixed_population_.data());

            // update the reference point
            //UpdateReferencePoint(g_GlobalSettings->mixed_population_.data(),2 * g_GlobalSettings->population_num_, g_GlobalSettings->obj_num_, reference_point_);

            // update ideal point
            UpdateIdealPoint(g_GlobalSettings->offspring_population_.data(),
                             g_GlobalSettings->population_num_, g_GlobalSettings->obj_num_, ideal_point_);

            // select next generation's population
            EnvironmentalSelection(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->mixed_population_.data());

            // record the population every interval generations and the first and last genration
            if (((g_GlobalSettings->output_interval_)&&
            (g_GlobalSettings->iteration_num_ % g_GlobalSettings->output_interval_ == 0))
            || g_GlobalSettings->IsTermination())
            {
                TrackPopulation(g_GlobalSettings->iteration_num_);
            }

            /**Py_BEGIN_ALLOW_THREADS, Py_END_ALLOW_THREADS
             * 为了在较长时间的C函数调用前，
             * 临时释放全局锁，
             * 完成后重新获取全局锁，
             * 以避免阻塞其他python的线程继续运行
             */
            // Py_END_ALLOW_THREADS;
        }
    }

    void R2_IBEA_LTR::Initialization()
    {
        // initialize parent population
        g_GlobalSettings->InitializePopulation(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);
        EvaluatePop(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);

        // set weight
        weight_ = SetWeight(g_GlobalSettings->weight_stringType_);
        // set step size
        step_size_ = g_GlobalSettings->step_size_;

        // load uniform weight vectors from *.txt
        char file[256];
        sprintf(file, "./UniformWeights/%dd_%d.txt",g_GlobalSettings->obj_num_, g_GlobalSettings->population_num_);
        lambda_ = LoadUniformWeights(g_GlobalSettings->population_num_, g_GlobalSettings->obj_num_, file);

        // initialize the reference point z* to be the utopian point (0, 0, ..., 0)
        InitializeReferencePoint(reference_point_, g_GlobalSettings->obj_num_);

        // initialize the ideal point
        UpdateIdealpoint(g_GlobalSettings->parent_population_.data(),
                         g_GlobalSettings->population_num_, ideal_point_, g_GlobalSettings->obj_num_);
    }

    void R2_IBEA_LTR::Crossover(Individual **parent_pop, Individual **offspring_pop)
    {
        // generate random permutation index for tournment selection
        int *index1 = new int[g_GlobalSettings->population_num_];
        int *index2 = new int[g_GlobalSettings->population_num_];
        random_permutation(index1, g_GlobalSettings->population_num_);
        random_permutation(index2, g_GlobalSettings->population_num_);

        for (int i = 0; i < g_GlobalSettings->population_num_ / 2; ++i)
        {
            Individual *parent1 = TournamentByBiR2Indicator(parent_pop[index1[2 * i]], parent_pop[index1[2 * i + 1]]);
            Individual *parent2 = TournamentByBiR2Indicator(parent_pop[index2[2 * i]], parent_pop[index2[2 * i + 1]]);
            SBX(parent1, parent2, offspring_pop[2 * i], offspring_pop[2 * i + 1], g_GlobalSettings);
        }
        delete[] index1;
        index1 = nullptr;
        delete[] index2;
        index2 = nullptr;
    }

    double R2_IBEA_LTR::CalEpsIndicator(Individual *ind1, Individual *ind2)
    {
        double max_eps = 0.0, temp_eps = 0.0;
        double r = g_GlobalSettings->dec_upper_bound_[0] - g_GlobalSettings->dec_lower_bound_[0];

        max_eps = (ind1->obj_[0] - g_GlobalSettings->dec_upper_bound_[0]) / r - (ind2->obj_[0] - g_GlobalSettings->dec_upper_bound_[0]) / r;
        for (int i = 1; i < g_GlobalSettings->obj_num_; ++i)
        {
            r = g_GlobalSettings->dec_upper_bound_[i] - g_GlobalSettings->dec_lower_bound_[1];
            temp_eps = (ind1->obj_[i] - g_GlobalSettings->dec_upper_bound_[i]) / r - (ind2->obj_[i] - g_GlobalSettings->dec_upper_bound_[i]) / r;

            if (temp_eps > max_eps)
                max_eps = temp_eps;
        }

        return max_eps;
    }

    void R2_IBEA_LTR::CalFitness(Individual **pop, int pop_num, double *fitness)
    {
        // determine indicator values and their maximum
        for (int i = 0; i < pop_num; ++i)
        {
            for (int j = 0; j < pop_num; ++j)
            {
                fitness[i * pop_num + j] = CalBiR2Indicator(pop[i], pop[j]);
            }
        }

        // calculate for each pair of individuals the corresponding value component
        for (int i = 0; i < pop_num;++i)
        {
            for (int j = 0; j < pop_num; j++)
            {
                fitness[i * pop_num + j] = -1 * exp(-1 * fitness[i * pop_num + j] / kappa_);
            }
        }

        // set individual's fitness
        for (int i = 0; i < pop_num; ++i)
        {
            double sum = 0;
            for (int j = 0; j < pop_num; ++j)
            {
                if (i != j)
                {
                    sum += fitness[j * pop_num + i];
                }
            }
            pop[i]->fitness_ = sum;
        }
    }

    void R2_IBEA_LTR::EnvironmentalSelection(Individual **parent_pop, Individual **mixed_pop)
    {
        int mixed_popnum = 2 * g_GlobalSettings->population_num_;

        // calculate fitness and store it in fitness array
        double *fitness = new double[mixed_popnum * mixed_popnum];
        CalFitness(mixed_pop, mixed_popnum, fitness);

        // select next generation's individuals
        int worst = -1;
        std::vector<int> flag(mixed_popnum,0);

        for (int i = g_GlobalSettings->population_num_; i > 0; --i)
        {
            int j = 0;
            for (j = 0; j < mixed_popnum && flag[j] == 1; ++j);
            worst = j;

            for (j = j + 1; j < mixed_popnum; ++j)
            {
                if (flag[j] != 1)
                {
                    if (mixed_pop[j]->fitness_ < mixed_pop[worst]->fitness_)
                        worst = j;
                }
            }

            for (j = 0; j < mixed_popnum; ++j)
            {
                if (flag[j] != 1)
                {
                    mixed_pop[j]->fitness_ -= fitness[worst * mixed_popnum + j];
                }
            }

            flag[worst] = 1;
        }

        int count = 0;
        for (int i = 0; i < mixed_popnum; ++i)
        {
            if (flag[i] != 1)
            {
                CopyIndividual(mixed_pop[i], parent_pop[count++]);
            }
        }
        if (fitness)
        {
            delete[] fitness;
            fitness= nullptr;
        }

    }

    double R2_IBEA_LTR::CalR2Indicator(Individual *x)
    {
        double res = 0;
        double a, ChebycheffValue;
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            ChebycheffValue = -INF;
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                a = fabs(ideal_point_[j] - x->obj_[j]) / lambda_[i][j];
                if (a > ChebycheffValue)
                {
                    ChebycheffValue = a;
                }
            }
            res += ChebycheffValue;
        }
        res /= g_GlobalSettings->population_num_;
        return res;
    }

    double R2_IBEA_LTR::CalR2Indicator(Individual *x, Individual *y)
    {
        double res = 0;
        double a, b, ChebycheffValue_x, ChebycheffValue_y;
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            ChebycheffValue_x = -INF;
            ChebycheffValue_y = -INF;
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                a = fabs(ideal_point_[j] - x->obj_[j]) / lambda_[i][j];
                b = fabs(ideal_point_[j] - y->obj_[j]) / lambda_[i][j];
                if (a > ChebycheffValue_x)
                {
                    ChebycheffValue_x = a;
                }
                if (b > ChebycheffValue_y)
                {
                    ChebycheffValue_y = b;
                }
            }
            if (ChebycheffValue_x < ChebycheffValue_y)
            {
                res += ChebycheffValue_x;
            }
            else
            {
                res += ChebycheffValue_y;
            }
        }
        res /= g_GlobalSettings->population_num_;
        return res;
    }

    double R2_IBEA_LTR::CalBiR2Indicator(Individual *x, Individual *y)
    {
        return CalR2Indicator(x) - CalR2Indicator(x, y);
    }

    Individual* R2_IBEA_LTR::TournamentByBiR2Indicator(Individual *ind1, Individual *ind2)
    {//TODO: debug
        double I_R2_1,I_R2_2;
        I_R2_1 = CalBiR2Indicator(ind1, ind2);
        I_R2_2 = CalBiR2Indicator(ind2, ind1);
        /*if (I_R2_1 < I_R2_2)
        {
            #if DEBUG
            std::cout<<"[INFO] Select parent 1."<<std::endl;
            #endif
            return ind1;
        }
        else if(I_R2_1 > I_R2_2)
        {
            #if DEBUG
            std::cout<<"[INFO] Select parent 2."<<std::endl;
            #endif
            return ind2;
        }*/
        if (I_R2_1 == 0)
        {
            return ind1;
        }
        else if(I_R2_2 == 0)
        {
            return ind2;
        }
        else
        {
            #if DEBUG
            std::cout<<"[INFO] Select parent randomly."<<std::endl;
            #endif
            if (randomperc() <= 0.5)
            {
                return (ind1);
            }
            else
            {
                return (ind2);
            }
        }
    }

    void R2_IBEA_LTR::UpdateTrainingSet()
    {
        double *temp, *ind_1_obj, *ind_2_obj;
        int index_1, index_2;
        int number_of_inquiries = g_GlobalSettings->inquiriesNum_;
        while(number_of_inquiries--)
        {
            if ((MAX_QUEUE_SIZE != -1) && (winners_.size()==MAX_QUEUE_SIZE))
            {
                temp = winners_.front();
                if (temp)
                {
                    delete[] temp;
                    temp = nullptr;
                }
                winners_.pop();
                temp = losers_.front();
                if (temp)
                {
                    delete[] temp;
                    temp = nullptr;
                }
                losers_.pop();
            }

            index_1 = index_2 = 0;
            while (index_1==index_2)
            {
                index_1 = rnd(0, g_GlobalSettings->population_num_ - 1);
                index_2 = rnd(0, g_GlobalSettings->population_num_ - 1);
            }
            ind_1_obj = new double[g_GlobalSettings->obj_num_];
            ind_2_obj = new double[g_GlobalSettings->obj_num_];
            for (int i = 0; i < g_GlobalSettings->obj_num_; ++i)
            {
                ind_1_obj[i] = g_GlobalSettings->parent_population_[index_1]->obj_[i];
                ind_2_obj[i] = g_GlobalSettings->parent_population_[index_2]->obj_[i];
            }
            if (CalInverseChebycheff(g_GlobalSettings->parent_population_[index_1], weight_, g_GlobalSettings->obj_num_)
            <CalInverseChebycheff(g_GlobalSettings->parent_population_[index_2], weight_, g_GlobalSettings->obj_num_))
            {
                winners_.push(ind_1_obj);
                losers_.push(ind_2_obj);
            }
            else
            {
                winners_.push(ind_2_obj);
                losers_.push(ind_1_obj);
            }
        }
    }

    void R2_IBEA_LTR::RecordCurrentPop(PyObject * pop)
    {
        PyObject *list;
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            list = PyList_New(g_GlobalSettings->obj_num_);//ref_cnt:=1
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                PyList_SetItem(list, j, PyFloat_FromDouble(g_GlobalSettings->parent_population_[i]->obj_[j]));
            }
            PyList_Append(pop, list);//ref_cnt:=2
            Py_DECREF(list);
        }
        // free space
        // PyList_SetSlice(list, 0, PyList_Size(list), NULL);
        // PySequence_DelSlice(list, 0, PySequence_Length(list));
    }

    void R2_IBEA_LTR::LoadTrainingSet(PyObject * winners, PyObject * losers)
    {
        PyObject *ind_obj_0, *ind_obj_1;
        for (int i = 0; i < winners_.size(); ++i)
        {
            ind_obj_0 = PyList_New(g_GlobalSettings->obj_num_);//ref_cnt:=1
            ind_obj_1 = PyList_New(g_GlobalSettings->obj_num_);//ref_cnt:=1
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                PyList_SetItem(ind_obj_0, j, PyFloat_FromDouble((winners_.front())[j]));
                PyList_SetItem(ind_obj_1, j, PyFloat_FromDouble((losers_.front())[j]));
            }
            winners_.push(winners_.front());
            winners_.pop();
            losers_.push(losers_.front());
            losers_.pop();
            PyList_Append(winners, ind_obj_0);//ref_cnt:=2
            PyList_Append(losers, ind_obj_1);//ref_cnt:=2
            Py_DECREF(ind_obj_0);
            Py_DECREF(ind_obj_1);
        }
        // free space
//        PyList_SetSlice(ind_obj_0, 0, PyList_Size(ind_obj_0), NULL);
//        PyList_SetSlice(ind_obj_1, 0, PyList_Size(ind_obj_1), NULL);
    }

    PyObject*  R2_IBEA_LTR::TrainRankNet_ReturnScore(PyObject *pFunction,PyObject *winners, PyObject *losers, PyObject *currPop)
    {
        PyObject *res;
        if (!pFunction)
        {
            std::cerr<<"[ERROR] pyFunction is NULL\n";
        }
        res = PyObject_CallFunctionObjArgs(pFunction, winners, losers, currPop, NULL);//res->ob_refcnt:=1
        std::cout<<"[INFO] PyObject_CallFunctionObjArgs(pFunc, winners_list, losers_list, PF_py, NULL) is over.\n";

        if (!res)
        {
            std::cerr<<"[ERROR] res is NULL, function in python excuted unsuccessfully."<<std::endl;
            return NULL;
        }
        else
        {
            std::cout<<"[INFO] Function excuted successfully."<<std::endl;
            return res;
        }
    }

    void R2_IBEA_LTR::UpdateScoreByRankNet(PyObject * res, double *score)
    {
        PyObject *item;
        int len = PyList_Size(res);
        for (int i = 0; i < len; ++i)
        {
            item = PyList_GetItem(res, i);
            score[i] = PyFloat_AsDouble(item);//item->ob_refcnt:=1
        }
    }

    void R2_IBEA_LTR::BiasingWeightSet()
    {
        #if UsingGoldenFunction
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            score_RankNet_[i] = 1 - CalInverseChebycheff(g_GlobalSettings->parent_population_[i], weight_, g_GlobalSettings->obj_num_);
        }
        #endif

        //score: the bigger, the better
        std::vector<double> scores(score_RankNet_,score_RankNet_+g_GlobalSettings->population_num_);
        std::vector<size_t> idx(scores.size());
        iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(),
             [&scores](size_t index_1, size_t index_2) { return scores[index_1] > scores[index_2]; });

        int solvedNum, currTuned, idx_of_nearest, currIndIdx, maxAttractionNum, nPromisingWeight;
        nPromisingWeight = (int)(g_GlobalSettings->population_num_*retention_rate_);
        maxAttractionNum = ceil((double)(g_GlobalSettings->population_num_- nPromisingWeight) / nPromisingWeight);
        int currPromisingWeightIndex[nPromisingWeight];
        double minDis, dis, minChebycheff, tempChebycheff;
        bool flag[g_GlobalSettings->population_num_];
        memset(flag, true, g_GlobalSettings->population_num_ * sizeof(bool));// true: adjustable, false: has been adjusted

        for (int i = 0; i < nPromisingWeight; ++i)
        {
            currIndIdx = idx[i];
            minChebycheff = INF;
            for (int j = 0; j < g_GlobalSettings->population_num_; ++j)
            {
                if (flag[j])
                {
                    tempChebycheff = CalInverseChebycheff(g_GlobalSettings->parent_population_[currIndIdx],
                                                          lambda_[j], ideal_point_, g_GlobalSettings->obj_num_);
                    if (tempChebycheff < minChebycheff)
                    {
                        minChebycheff = tempChebycheff;
                        currPromisingWeightIndex[i] = j;
                    }
                }
            }
            flag[currPromisingWeightIndex[i]] = false;
        }
        solvedNum = nPromisingWeight;

        for (int i = 0; (i < nPromisingWeight) && (solvedNum < g_GlobalSettings->population_num_); ++i)
        {
            currTuned = 0;
            while (currTuned<maxAttractionNum && (solvedNum < g_GlobalSettings->population_num_))
            {
                minDis = INF;
                //find the nearest weight vector from the remaining set of weight vectors
                for (int j = 0; j < g_GlobalSettings->population_num_; ++j)
                {
                    if (flag[j])
                    {
                        dis = CalEuclidianDistance(lambda_[currPromisingWeightIndex[i]], lambda_[j], g_GlobalSettings->obj_num_);
                        if (dis < minDis)
                        {
                            minDis = dis;
                            idx_of_nearest = j;
                        }
                    }
                }
                for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
                {
                    lambda_[idx_of_nearest][j] += step_size_ * (lambda_[currPromisingWeightIndex[i]][j] - lambda_[idx_of_nearest][j]);
                }
                flag[idx_of_nearest] = false;//solved
                solvedNum++;
                currTuned++;
            }
        }
        #if DEBUG
        std::cout<<"[INFO] Adjusted weight vectors: \n";
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                std::cout<<lambda_[i][j]<<"\t";
            }
            std::cout<<std::endl;
        }
        #endif
    }

    void R2_IBEA_LTR::Consultation_PreferenceElicitation()
    {
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append('./')");
        wchar_t * s2[] = { L" " }; // 宽字符，长度为2字节
        PySys_SetArgv(1, s2);   // TODO: 加入argv参数 否则出错.PyAPI_FUNC(void) PySys_SetArgv(int, wchar_t **);

        PyObject *pModule = NULL;
        PyObject *pTrainModel = NULL;
        PyObject *pUseModel = NULL;
        PyObject *pRet = NULL;
        PyObject *winners_list = NULL;
        PyObject *losers_list = NULL;
        PyObject *currPop_list = NULL;

        //导入python文件
        pModule = PyImport_ImportModule("RankNet");
        if (!pModule)
        {
            std::cerr << "[ERROR] Can not open python file!\n";
            return;
        }

        /***
         * initialize model in first consultation
         */
        if (g_GlobalSettings->iteration_num_ == g_GlobalSettings->tau_)
        {
            PyObject *pCreateModel = PyObject_GetAttrString(pModule, "CreateModel");
            PyObject_CallFunctionObjArgs(pCreateModel, PyLong_FromLong((long)(g_GlobalSettings->obj_num_)),PyLong_FromLong((long)(g_GlobalSettings->run_id_)), NULL);
            Py_DECREF(pCreateModel);
        }

        winners_list = PyList_New(0);
        losers_list = PyList_New(0);
        currPop_list = PyList_New(0);

        RecordCurrentPop(currPop_list);
        UpdateTrainingSet();
        LoadTrainingSet(winners_list, losers_list);

        pTrainModel = PyObject_GetAttrString(pModule, "TrainModel");
        if (!pTrainModel || !PyCallable_Check(pTrainModel))
        {
            std::cerr << "[ERROR] Can't find function" << std::endl;
            return;
        }
        PyObject_CallFunctionObjArgs(pTrainModel, winners_list,losers_list,PyLong_FromLong((long)(g_GlobalSettings->run_id_)), NULL);
        pUseModel = PyObject_GetAttrString(pModule, "UseModel");
        pRet = PyObject_CallFunctionObjArgs(pUseModel, currPop_list,PyLong_FromLong((long)(g_GlobalSettings->run_id_)), NULL);
        // in case failed
        if (!pRet)
        {
            std::cerr<<"[ERROR] TrainRankNet_ReturnScore failed.\n";
            return;
        }
        else
        {
            // assign scores to current population via trained RankNet
            UpdateScoreByRankNet(pRet, score_RankNet_);
        }
        BiasingWeightSet();

        //int nlist=PyList_ClearFreeList();
        /*PyList_SetSlice(winners_list, 0, PyList_Size(winners_list), NULL);
        PyList_SetSlice(losers_list, 0, PyList_Size(losers_list), NULL);
        PyList_SetSlice(currPop_list, 0, PyList_Size(currPop_list), NULL);*/

        #if DEBUG
        //std::cout<<"[INFO] Return of PyList_ClearFreeList: "<<nlist<<std::endl;
        std::cout<<"[INFO] pRet->ob_refcnt: "<<pRet->ob_refcnt<<std::endl;
        std::cout<<"[INFO] currPop_list->ob_refcnt: "<<currPop_list->ob_refcnt<<std::endl;
        std::cout<<"[INFO] winners_list->ob_refcnt: "<<winners_list->ob_refcnt<<std::endl;
        std::cout<<"[INFO] losers_list->ob_refcnt: "<<losers_list->ob_refcnt<<std::endl;
        std::cout<<"[INFO] pModule->ob_refcnt: "<<pModule->ob_refcnt<<std::endl;
        #endif

        /**
         * free list space
         */
        /*PyObject *item;
        int len = PyList_Size(currPop_list);
        std::cout<<"len:="<<len<<std::endl;
        for (int i = 0; i < len; ++i)
        {
            item = PyList_GetItem(currPop_list, i);
            std::cout<<"[INFO] item->ob_refcnt: "<<item->ob_refcnt<<std::endl;
            Py_DECREF(item);
        }
        len = PyList_Size(winners_list);
        std::cout<<"len:="<<len<<std::endl;
        for (int i = 0; i < len; ++i)
        {
            item = PyList_GetItem(winners_list, i);
            std::cout<<"[INFO] item->ob_refcnt: "<<item->ob_refcnt<<std::endl;
            Py_DECREF(item);
            item = PyList_GetItem(losers_list, i);
            std::cout<<"[INFO] item->ob_refcnt: "<<item->ob_refcnt<<std::endl;
            Py_DECREF(item);
        }*/

        #if DEBUG
        //std::cout<<"[INFO] Return of PyList_ClearFreeList: "<<nlist<<std::endl;
        std::cout<<"[INFO] pRet->ob_refcnt: "<<pRet->ob_refcnt<<std::endl;
        std::cout<<"[INFO] currPop_list->ob_refcnt: "<<currPop_list->ob_refcnt<<std::endl;
        std::cout<<"[INFO] winners_list->ob_refcnt: "<<winners_list->ob_refcnt<<std::endl;
        std::cout<<"[INFO] losers_list->ob_refcnt: "<<losers_list->ob_refcnt<<std::endl;
        std::cout<<"[INFO] pModule->ob_refcnt: "<<pModule->ob_refcnt<<std::endl;
        #endif

        Py_DECREF(pModule);
        Py_DECREF(pTrainModel);//TODO: DEBUG -> ref_cnt != 0
        Py_DECREF(pUseModel);
        Py_DECREF(pRet);
        Py_DECREF(currPop_list);
        Py_DECREF(winners_list);
        Py_DECREF(losers_list);
    }

}
