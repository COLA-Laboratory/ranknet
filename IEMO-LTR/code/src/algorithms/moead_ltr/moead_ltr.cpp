//
// Created by guiyulai on 2021/6/9.
//

#include "algorithms/moead_ltr/moead_ltr.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <numeric>
#include <cstring>

#include "core/global.h"
#include "core/utility.h"
#include "core/uniform_point.h"
#include "operator/sbx.h"
#include "operator/mutation.h"
#include "random/random.h"
using namespace std;



namespace emoc {
    static int all_comparisons=0;
    static int errors=0;

    MOEAD_LTR::MOEAD_LTR(Problem *problem, int thread_num):
            Algorithm(problem,thread_num),
            lambda_(nullptr),
            weight_(nullptr),
            goldenPoint_(nullptr),
            weight_num_(0),
            neighbour_(nullptr),
            ideal_point_(new double[g_GlobalSettings->obj_num_]),
            aggregation_type_(0),
            pbi_theta_(5.0),
            step_size_(0.5),
            retention_rate_(0.2),
            score_RankNet_(new double[g_GlobalSettings->population_num_])
    {

    }

    MOEAD_LTR::~MOEAD_LTR()
    {
        for (int i = 0; i < weight_num_; ++i)
        {
            delete[] lambda_[i];
            delete[] neighbour_[i];
            lambda_[i] = nullptr;
            neighbour_[i] = nullptr;
        }
        delete[] lambda_;
        delete[] goldenPoint_;
        delete[] neighbour_;
        delete[] ideal_point_;
        delete[] weight_;
        delete[] score_RankNet_;
        lambda_ = nullptr;
        goldenPoint_ = nullptr;
        neighbour_ = nullptr;
        ideal_point_ = nullptr;
        weight_ = nullptr;
        score_RankNet_ = nullptr;
        while(winners_.size())
        {
            delete[] winners_.front();
            winners_.pop();
            delete[] losers_.front();
            losers_.pop();
        }
    }

    void MOEAD_LTR::Run()
    {
        // TODO: multi-thread for Python
        class PyThreadStateLock PyThreadLock;

        Initialization();
        #ifdef GNUPLOT_COMMAND
        FILE *gp= nullptr;
        gp = popen(GNUPLOT_COMMAND,"w");
        if (!gp)
        {
            cerr<<"<Error!!!> Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND"<<endl;
            exit(-1);
        }
        #endif
        if (g_GlobalSettings->output_interval_)
        {
            TrackPopulation(g_GlobalSettings->iteration_num_);
        }
        Individual *offspring = g_GlobalSettings->offspring_population_[0];

        while (!g_GlobalSettings->IsTermination())
        {
            #ifdef GNUPLOT_COMMAND
            display_pop(gp,g_GlobalSettings->parent_population_.data(),
                        g_GlobalSettings->population_num_,g_GlobalSettings->obj_num_,
                        g_GlobalSettings->iteration_num_);
            #endif
            // begin each iteration
            g_GlobalSettings->iteration_num_++;
            // consultation, preference elicitation(tune the weights)
            if (g_GlobalSettings->tau_&&(g_GlobalSettings->iteration_num_ % g_GlobalSettings->tau_ == 0))
            {// if g_GlobalSettings->tau_ is equal to 0, means no consultation and preference elicitation
                #if UsingGoldenFunction
                UsingGoldenFunc();
                #else
                Consultation_PreferenceElicitation();
                #endif
            }

            /**Py_BEGIN_ALLOW_THREADS, Py_END_ALLOW_THREADS
             * 为了在较长时间的C函数调用前，
             * 临时释放全局锁，
             * 完成后重新获取全局锁，
             * 以避免阻塞其他python的线程继续运行
             */
             Py_BEGIN_ALLOW_THREADS;

            UpdateNeighbours();

            for (int i = 0; i < weight_num_; ++i)
            {
                // generate offspring for current subproblem
                Crossover(g_GlobalSettings->parent_population_.data(), i, offspring);
                MutationInd(offspring, g_GlobalSettings);
                EvaluateInd(offspring);//g_GlobalSettings->current_evaluation_++;

                // update ideal point
                UpdateIdealpoint(offspring, ideal_point_, g_GlobalSettings->obj_num_);

                // update neighbours' subproblem
                UpdateSubproblem(offspring, i, aggregation_type_);
            }

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
             Py_END_ALLOW_THREADS;
        }
        #ifdef GNUPLOT_COMMAND
        pclose(gp);
        #endif
    }

    void MOEAD_LTR::Initialization()
    {
        weight_num_ = g_GlobalSettings->population_num_;
        // initialize parent population
        g_GlobalSettings->InitializePopulation(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);
        EvaluatePop(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);

        // set weight
        weight_ = SetWeight(g_GlobalSettings->weight_stringType_);
        step_size_ = g_GlobalSettings->step_size_;

        // load uniform weights from *.txt generated by hitandrun
        char file[256];
        sprintf(file, "./UniformWeights/%dd_%d.txt",g_GlobalSettings->obj_num_, g_GlobalSettings->population_num_);
        lambda_ = LoadUniformWeights(g_GlobalSettings->population_num_, g_GlobalSettings->obj_num_, file);

        // set the neighbours of each individual
        InitializeNeighbours();

        // initialize ideal point
        UpdateIdealpoint(g_GlobalSettings->parent_population_.data(), weight_num_, ideal_point_, g_GlobalSettings->obj_num_);
    }

    void MOEAD_LTR::InitializeNeighbours()
    {
        // set neighbour size and allocate memory
        neighbour_num_ = weight_num_ / 10;//todo: neighbour size
        neighbour_ = new int*[weight_num_];
        for (int i = 0; i < weight_num_; ++i)
        {
            neighbour_[i] = new int[neighbour_num_];
        }

        DistanceInfo *sort_list = new DistanceInfo[weight_num_];
        for (int i = 0; i < weight_num_; ++i)
        {
            for (int j = 0; j < weight_num_; ++j)
            {
                // calculate distance to each weight vector
                double distance_temp = 0;
                for (int k = 0; k < g_GlobalSettings->obj_num_; ++k)
                {
                    distance_temp += (lambda_[i][k] - lambda_[j][k]) * (lambda_[i][k] - lambda_[j][k]);
                }

                sort_list[j].distance = sqrt(distance_temp);
                sort_list[j].index = j;
            }

            std::sort(sort_list, sort_list+weight_num_, [](const DistanceInfo &left, const DistanceInfo &right) {
                return left.distance < right.distance;
            });

            for (int j = 0; j < neighbour_num_; j++)
            {
                neighbour_[i][j] = sort_list[j+1].index;
            }
        }

        delete[] sort_list;
    }

    void MOEAD_LTR::UpdateNeighbours()
    {
        DistanceInfo *sort_list = new DistanceInfo[weight_num_];
        for (int i = 0; i < weight_num_; ++i)
        {
            for (int j = 0; j < weight_num_; ++j)
            {
                // calculate distance to each weight vector
                double distance_temp = 0;
                for (int k = 0; k < g_GlobalSettings->obj_num_; ++k)
                {
                    distance_temp += (lambda_[i][k] - lambda_[j][k]) * (lambda_[i][k] - lambda_[j][k]);
                }

                sort_list[j].distance = sqrt(distance_temp);
                sort_list[j].index = j;
            }

            std::sort(sort_list, sort_list+weight_num_, [](const DistanceInfo &left, const DistanceInfo &right) {
                return left.distance < right.distance;
            });

            for (int j = 0; j < neighbour_num_; j++)
            {
                neighbour_[i][j] = sort_list[j+1].index;
            }
        }
        delete[] sort_list;
    }

    void MOEAD_LTR::Crossover(Individual **parent_pop, int current_index, Individual *offspring)
    {
        // randomly select two parent from current individual's neighbours
        int k = rnd(0, neighbour_num_ - 1);
        int l = rnd(0, neighbour_num_ - 1);
        Individual *parent1 = parent_pop[neighbour_[current_index][k]];
        Individual *parent2 = parent_pop[neighbour_[current_index][l]];
        SBX(parent1, parent2, g_GlobalSettings->offspring_population_[1], offspring, g_GlobalSettings);
    }

    void MOEAD_LTR::UpdateSubproblem(Individual *offspring, int current_index, int aggregation_type)
    {
        double *offspring_fitness = new double[neighbour_num_];
        double *neighbour_fitness = new double[neighbour_num_];

        // calculate fitness;
        switch (aggregation_type)
        {
            case 0:
                // inverse chebycheff
                for (int i = 0; i < neighbour_num_; ++i)
                {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = CalInverseChebycheff(current_ind, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_);
                    offspring_fitness[i] = CalInverseChebycheff(offspring, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_);
                }
                break;

            case 1:
                // weighted sum
                for (int i = 0; i < neighbour_num_; ++i)
                {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = CalWeightedSum(current_ind, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_);
                    offspring_fitness[i] = CalWeightedSum(offspring, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_);
                }
                break;

            case 2:
                // PBI
                for (int i = 0; i < neighbour_num_; ++i)
                {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = CalPBI(current_ind, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_, pbi_theta_);
                    offspring_fitness[i] = CalPBI(offspring, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_, pbi_theta_);
                }
                break;

            default:
                break;
        }

        // update subproblem
        for (int i = 0; i < neighbour_num_; ++i)
        {
            if (offspring_fitness[i] < neighbour_fitness[i])
            {
                CopyIndividual(offspring, g_GlobalSettings->parent_population_[neighbour_[current_index][i]]);
            }
        }


        delete[] neighbour_fitness;
        delete[] offspring_fitness;
    }
    void MOEAD_LTR:: Consultation_PreferenceElicitation()
    {
        PyObject *pRet = NULL;
        PyObject *winners_list = NULL;
        PyObject *losers_list = NULL;
        PyObject *currPop_list = NULL;
        PyObject *pModule = NULL;
        PyObject *pTrainModel = NULL;
        PyObject *pUseModel = NULL;

        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append('./')");
        wchar_t * s2[] = { L" " }; // 宽字符，长度为2字节
        PySys_SetArgv(1, s2);   // TODO: 加入argv参数 否则出错.PyAPI_FUNC(void) PySys_SetArgv(int, wchar_t **);
        //导入python文件
        pModule = PyImport_ImportModule("RankNet");
        if (!pModule)
        {
            std::cerr << "[ERROR] Can not open python file!\n";
            return;
        }
        DEBUG_MSG("pModule->ob_refcnt at the beginning",pModule->ob_refcnt);


        /***
         * initialize model in first consultation
         */
        if (g_GlobalSettings->iteration_num_ == g_GlobalSettings->tau_)
        {
            PyObject *pCreateModel = PyObject_GetAttrString(pModule, "CreateModel");
            DEBUG_MSG("pCreateModel->ob_refcnt before Call",pCreateModel->ob_refcnt);
            PyObject_CallFunctionObjArgs(pCreateModel, PyLong_FromLong((long)(g_GlobalSettings->obj_num_)),PyLong_FromLong((long)(g_GlobalSettings->run_id_)), NULL);
            DEBUG_MSG("pCreateModel->ob_refcnt",pCreateModel->ob_refcnt);
            Py_DECREF(pCreateModel);
            DEBUG_MSG("pCreateModel->ob_refcnt",pCreateModel->ob_refcnt);
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
        PyObject_CallFunctionObjArgs(pTrainModel,winners_list,losers_list,PyLong_FromLong((long)(g_GlobalSettings->run_id_)), NULL);
        pUseModel = PyObject_GetAttrString(pModule, "UseModel");
        pRet = PyObject_CallFunctionObjArgs(pUseModel, currPop_list,PyLong_FromLong((long)(g_GlobalSettings->run_id_)), NULL);
//        pRet = TrainRankNet_ReturnScore(pFunc,winners_list, losers_list, currPop_list);
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

        /**
         * free list space
         */
        DEBUG_MSG("pModule->ob_refcnt",pModule->ob_refcnt);
        DEBUG_MSG("pTrainModel->ob_refcnt",pTrainModel->ob_refcnt);
        DEBUG_MSG("pUseModel->ob_refcnt",pUseModel->ob_refcnt);
        DEBUG_MSG("pRet->ob_refcnt",pRet->ob_refcnt);
        DEBUG_MSG("currPop_list->ob_refcnt",currPop_list->ob_refcnt);
        DEBUG_MSG("winners_list->ob_refcnt",winners_list->ob_refcnt);
        DEBUG_MSG("losers_list->ob_refcnt",losers_list->ob_refcnt);

        Py_DECREF(pModule);
        Py_DECREF(pTrainModel);
        Py_DECREF(pUseModel);
        Py_DECREF(pRet);
        Py_DECREF(currPop_list);
        Py_DECREF(winners_list);
        Py_DECREF(losers_list);
    }

    void MOEAD_LTR::UsingGoldenFunc()
    {
        std::cout<<"<info> using golden function as preference"<<endl;
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            score_RankNet_[i] = 1 - CalInverseChebycheff(g_GlobalSettings->parent_population_[i], weight_, g_GlobalSettings->obj_num_);
        }
        #if DEBUG
        // track ranking performance of top-k individuals
        double goldenPreference[g_GlobalSettings->population_num_];
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            goldenPreference[i] = CalInverseChebycheff(g_GlobalSettings->parent_population_[i]->obj_, weight_, g_GlobalSettings->obj_num_);
        }
        double acc=topK(score_RankNet_, goldenPreference, g_GlobalSettings->population_num_, (int)(g_GlobalSettings->population_num_*0.2));
        DEBUG_MSG("Accuracy of top-k ranking", acc);
        #endif
        BiasingWeightSet();
    }
    void MOEAD_LTR::RecordCurrentPop(PyObject *pop)
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
    }
    void MOEAD_LTR::UpdateTrainingSet()
    {
        double *temp;
        double *ind_1_obj;
        double *ind_2_obj;
        double golden_function[2];
        int index_1, index_2;
        int number_of_inquiries = g_GlobalSettings->inquiriesNum_;

        while(number_of_inquiries--)
        {
            ++all_comparisons;
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
            ind_1_obj=new double[g_GlobalSettings->obj_num_];
            ind_2_obj=new double[g_GlobalSettings->obj_num_];
            for (int i = 0; i < g_GlobalSettings->obj_num_; ++i)
            {
                ind_1_obj[i] = g_GlobalSettings->parent_population_[index_1]->obj_[i];
                ind_2_obj[i] = g_GlobalSettings->parent_population_[index_2]->obj_[i];
            }
            golden_function[0]=CalInverseChebycheff(g_GlobalSettings->parent_population_[index_1], weight_, g_GlobalSettings->obj_num_);
            golden_function[1]=CalInverseChebycheff(g_GlobalSettings->parent_population_[index_2], weight_, g_GlobalSettings->obj_num_);

            if (g_GlobalSettings->kappa_<0)
            {
                if (golden_function[0]<golden_function[1])
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
            else//add some noise errors
            {
                double delta=fabs(golden_function[0]-golden_function[1]);
                DEBUG_MSG("delta",delta)
                double error_prob=exp(-1.0*g_GlobalSettings->kappa_*delta);
                DEBUG_MSG("error_prob",error_prob)
                if (randomperc()<error_prob)// randomly make error
                {
                    ++errors;
                    DEBUG_MSG("randomly make error","")
                    DEBUG_MSG("errors proportion so far",(double)(errors)/(double)(all_comparisons))
                    if (golden_function[0]<golden_function[1])
                    {
                        winners_.push(ind_2_obj);
                        losers_.push(ind_1_obj);
                    }
                    else
                    {
                        winners_.push(ind_1_obj);
                        losers_.push(ind_2_obj);
                    }

                }
                else// randomly make no error
                {
                    DEBUG_MSG("randomly make no error","")
                    DEBUG_MSG("no errors proportion so far",(double)(all_comparisons-errors)/(double)(all_comparisons))
                    if (golden_function[0]<golden_function[1])
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
        }
    }
    void MOEAD_LTR:: LoadTrainingSet(PyObject *winners, PyObject *losers)
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
    }
    PyObject*  MOEAD_LTR::TrainRankNet_ReturnScore(PyObject *pFunction,PyObject *winners, PyObject *losers, PyObject *currPop)
    {
        PyObject *res;
        if (!pFunction)
        {
            std::cerr<<"[ERROR] pyFunction is NULL\n";
        }
        res = PyObject_CallFunctionObjArgs(pFunction, winners, losers, currPop, NULL);//res->ob_refcnt:=1
        if (!res)
        {
            std::cerr<<"[ERROR] res is NULL, function in python excuted unsuccessfully."<<std::endl;
            return NULL;
        }
        else
        {
            return res;
        }
    }
    void MOEAD_LTR:: UpdateScoreByRankNet(PyObject * res, double *score)
    {
        PyObject *item;
        int len = PyList_Size(res);
        for (int i = 0; i < len; ++i)
        {
            item = PyList_GetItem(res, i);
            score[i] = PyFloat_AsDouble(item);//item->ob_refcnt:=1
        }
        #if DEBUG
        // track ranking performance of top-k individuals
        double goldenPreference[g_GlobalSettings->population_num_];
        for (int i = 0; i < g_GlobalSettings->population_num_; ++i)
        {
            goldenPreference[i] = CalInverseChebycheff(g_GlobalSettings->parent_population_[i]->obj_, weight_, g_GlobalSettings->obj_num_);
        }
        double acc=topK(score, goldenPreference, g_GlobalSettings->population_num_, (int)(g_GlobalSettings->population_num_*0.2));
        DEBUG_MSG("Accuracy of top-k ranking", acc);
        #endif
    }
    void MOEAD_LTR:: BiasingWeightSet()
    {//TODO: adjust weight vectors associated with promising weight
        //score: the bigger, the better
        std::vector<double> scores(score_RankNet_,score_RankNet_+g_GlobalSettings->population_num_);
        std::vector<size_t> idx(scores.size());
        iota(idx.begin(), idx.end(), 0);
        sort(idx.begin(), idx.end(),
             [&scores](size_t index_1, size_t index_2) { return scores[index_1] > scores[index_2]; });

        int solvedNum, currTuned, idx_of_nearest, currIndIdx, maxAttractionNum, nPromisingWeight;
        nPromisingWeight = (int)(g_GlobalSettings->population_num_*retention_rate_);
        maxAttractionNum = ceil((double)(g_GlobalSettings->population_num_- nPromisingWeight) / nPromisingWeight);
//        int currPromisingWeightIndex[nPromisingWeight];
        double minDis, dis, minChebycheff, tempChebycheff;
        bool flag[g_GlobalSettings->population_num_];
        memset(flag, true, g_GlobalSettings->population_num_ * sizeof(bool));// true: adjustable, false: has been adjusted

        // TODO: finding promising weight vectors of promising solutions using Tchebycheff rather than associated weight vectors
        /*for (int i = 0; i < nPromisingWeight; ++i)
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
        }*/
        for (int i = 0; i < nPromisingWeight; ++i)
        {
            flag[idx[i]] = false;
        }

        #if DEBUG
        /*std::cout<<"<info> currPromisingWeightIndex:\n";
        for (int i = 0; i < nPromisingWeight; ++i)
        {
            std::cout<<idx[i]<<"\t"<<currPromisingWeightIndex[i]<<std::endl;
        }*/
        #endif
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
                        dis = CalEuclidianDistance(lambda_[idx[i]], lambda_[j], g_GlobalSettings->obj_num_);
                        if (dis < minDis)
                        {
                            minDis = dis;
                            idx_of_nearest = j;
                        }
                    }
                }
                for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
                {
                    lambda_[idx_of_nearest][j] += step_size_ * (lambda_[idx[i]][j] - lambda_[idx_of_nearest][j]);
                }
                flag[idx_of_nearest] = false;//solved
                solvedNum++;
                currTuned++;
            }
        }
    }

}
