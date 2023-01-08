//
// Created by guiyulai on 2021/6/8.
//


#include "algorithms/plvf/plvf.h"

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

    PLVF::PLVF(Problem *problem, int thread_num) :
            Algorithm(problem, thread_num),
            lambda_(nullptr),
            weight_(nullptr),
            goldenPoint_(nullptr),
            weight_num_(0),
            neighbour_(nullptr),
            ideal_point_(new double[g_GlobalSettings->obj_num_]),
            aggregation_type_(0),
            pbi_theta_(5.0),
            alpha_(5.0),
            alpha_DM_(5.0),
            rbf_output_(nullptr),
            rbfnet_c_(nullptr),
            rbfnet_weight_(nullptr),
            mu_(0) ,
            sigma_(20),
            stepSize_(0.5)
            {}

    PLVF::~PLVF() {
        for (int i = 0; i < weight_num_; ++i)
        {
            delete[] lambda_[i];
            delete[] neighbour_[i];
            lambda_[i] = nullptr;
            neighbour_[i] = nullptr;
        }
        for (int i = 0; i < 40; ++i)
        {
            delete[] rbfnet_c_[i];
            rbfnet_c_[i] = nullptr;
        }
        delete[] lambda_;
        delete[] weight_;
        delete[] goldenPoint_;
        delete[] neighbour_;
        delete[] ideal_point_;
        delete[] rbf_output_;
        delete[] rbfnet_weight_;
        delete[] rbfnet_c_;
        lambda_ = nullptr;
        weight_ = nullptr;
        goldenPoint_ = nullptr;
        neighbour_ = nullptr;
        ideal_point_ = nullptr;
        rbf_output_ = nullptr;
        rbfnet_weight_ = nullptr;
        rbfnet_c_ = nullptr;
    }

    void PLVF::Run() {
        Initialization();
        Individual *offspring = g_GlobalSettings->offspring_population_[0];
        int maxGen = (int) (g_GlobalSettings->max_evaluation_ / g_GlobalSettings->population_num_);//TODO: DEBUG
        int tau;
        if (g_GlobalSettings->tau_)
        {
            tau = g_GlobalSettings->tau_;
        }
        else
        {
            tau = ceil((double) maxGen / (double) (g_GlobalSettings->askTime_ + 1));
        }

        int number_of_inquiries;
        int randNum;
        SortList *fitnessList = new SortList[weight_num_];
        double dis, minDis;
        int idx_of_nearest, currPromisingIndex, maxAttractionNum, currTuned;
        bool flag[weight_num_];
        int solvedNum;

        if (g_GlobalSettings->output_interval_)
        {
            TrackPopulation(g_GlobalSettings->iteration_num_);
        }
        while (!g_GlobalSettings->IsTermination()) {
            // begin each iteration
            g_GlobalSettings->iteration_num_++;
            // do consultation and tune the weights
            if (g_GlobalSettings->iteration_num_ % tau == 0)
            {
                if (g_GlobalSettings->iteration_num_ == tau)
                {//select randomly for the first consultation
                    number_of_inquiries = mu_;
                    for (int i = 0; i < number_of_inquiries; ++i)
                    {
                        randNum = rnd(0, weight_num_ - 1);
                        rbf_output_[i] = CalInverseChebycheff(g_GlobalSettings->parent_population_[randNum], weight_,
                                                              ideal_point_, g_GlobalSettings->obj_num_);
                        for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
                        {
                            rbfnet_c_[i][j] = g_GlobalSettings->parent_population_[randNum]->obj_[j];
                        }
                    }
                }
                else
                {//select the best top-10 individuals by RBFNet
                    number_of_inquiries = g_GlobalSettings->inquiriesNum_;
                    for (int i = 0; i < weight_num_; ++i)
                    {
                        UsingRBFNet(g_GlobalSettings->parent_population_[i], rbfnet_c_, rbfnet_weight_, sigma_, number_of_inquiries, g_GlobalSettings->obj_num_);
                        fitnessList[i].value = g_GlobalSettings->parent_population_[i]->fitness_;
                        fitnessList[i].index = i;
                    }
                    std::sort(fitnessList, fitnessList+weight_num_, [](const SortList &left, const SortList &right) {
                        return left.value < right.value;
                    });
                    solvedNum = number_of_inquiries;
                    memset(flag, true, weight_num_ * sizeof(bool));// true: adjustable; false: has been adjusted
                    for (int i = 0; i < number_of_inquiries; ++i)
                    {
                        flag[fitnessList[i].index] = false;
                        rbf_output_[i] = CalInverseChebycheff(g_GlobalSettings->parent_population_[fitnessList[i].index], weight_,
                                                              ideal_point_, g_GlobalSettings->obj_num_);
                        for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
                        {
                            rbfnet_c_[i][j] = g_GlobalSettings->parent_population_[fitnessList[i].index]->obj_[j];
                        }
                    }
                    maxAttractionNum = ceil((double)(weight_num_- number_of_inquiries) / number_of_inquiries);
                    for (int i = 0; (i<number_of_inquiries)&&(solvedNum<weight_num_); ++i)
                    {
                        currPromisingIndex = fitnessList[i].index;
                        currTuned = 0;
                        while (currTuned<maxAttractionNum && (solvedNum < weight_num_))
                        {
                            minDis = INF;
                            //find the nearest weight vector from the remaining set of weight vectors
                            for (int j = 0; j < weight_num_; ++j)
                            {
                                if (flag[j])
                                {
                                    dis = CalEuclidianDistance(lambda_[currPromisingIndex], lambda_[j], g_GlobalSettings->obj_num_);//TODO: consine?
                                    if (dis < minDis)
                                    {
                                        minDis = dis;
                                        idx_of_nearest = j;
                                    }
                                }
                            }
                            //cout<<"tune w["<<idx_of_nearest<<"] to w["<<currPromisingIndex<<"]"<<endl;
                            //move the closest weight vector to the current promising weight vector
                            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
                            {
                                lambda_[idx_of_nearest][j] += stepSize_ * (lambda_[currPromisingIndex][j] - lambda_[idx_of_nearest][j]);
                            }
                            flag[idx_of_nearest] = false;//solved
                            solvedNum++;
                            currTuned++;
                        }
                    }

                }
                TrainRBFNet(rbfnet_c_, rbfnet_weight_, rbf_output_, sigma_, number_of_inquiries,
                            g_GlobalSettings->obj_num_);
                // Update the neighborhood
                UpdateNeighbours();
            }

            for (int i = 0; i < weight_num_; ++i) {
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
        }
    }

    void PLVF::Initialization() {
        weight_num_ = g_GlobalSettings->population_num_;
        // TODO: rbfnet memory allocation
        rbf_output_ = new double[weight_num_];
        rbfnet_c_ = new double*[40];
        for (int i = 0; i < 40; i++) {
            rbfnet_c_[i] = new double[g_GlobalSettings->obj_num_];
        }
        rbfnet_weight_ = new double[40];
        mu_ = 2 * g_GlobalSettings->obj_num_ + 1;
        if (mu_ > 40)
        {
            cerr << "mu_>40" << endl;
            exit(1);
        }

        // initialize parent population
        g_GlobalSettings->InitializePopulation(g_GlobalSettings->parent_population_.data(),
                                               weight_num_);
        EvaluatePop(g_GlobalSettings->parent_population_.data(), weight_num_);

        // set weight
        weight_ = SetWeight(g_GlobalSettings->weight_stringType_);
        stepSize_ = g_GlobalSettings->step_size_;

        // load unifrom weights from *.txt generated by hitandrun
        char file[256];
        sprintf(file, "./UniformWeights/%dd_%d.txt",g_GlobalSettings->obj_num_, g_GlobalSettings->population_num_);
        lambda_ = LoadUniformWeights(g_GlobalSettings->population_num_, g_GlobalSettings->obj_num_, file);

        // set the neighbours of each individual
        InitializeNeighbours();

        // initialize ideal point
        UpdateIdealpoint(g_GlobalSettings->parent_population_.data(), weight_num_, ideal_point_,
                         g_GlobalSettings->obj_num_);

    }

    void PLVF::InitializeNeighbours() {
        // set neighbour size and allocate memory
        neighbour_num_ = 10;//todo: neighborhood size
        neighbour_ = new int *[weight_num_];
        for (int i = 0; i < weight_num_; ++i)
        {
            neighbour_[i] = new int[neighbour_num_];
        }
        DistanceInfo *sort_list = new DistanceInfo[weight_num_];
        for (int i = 0; i < weight_num_; ++i) {
            for (int j = 0; j < weight_num_; ++j) {
                /*// calculate distance to each weight vector
                double distance_temp = 0;
                for (int k = 0; k < g_GlobalSettings->obj_num_; ++k)
                {
                    distance_temp += (lambda_[i][k] - lambda_[j][k]) * (lambda_[i][k] - lambda_[j][k]);
                }
                sort_list[j].distance = sqrt(distance_temp);
                sort_list[j].index = j;*/
                /**
                 * find neighbors by cosine, not euclidean disdance
                 */
                sort_list[j].distance = CalculateCos(lambda_[i], lambda_[j], g_GlobalSettings->obj_num_);
                sort_list[j].index = j;

            }
            std::sort(sort_list, sort_list + weight_num_, [](const DistanceInfo &left, const DistanceInfo &right) {
                return left.distance > right.distance;
            });
            for (int j = 0; j < neighbour_num_; j++) {
                neighbour_[i][j] = sort_list[j + 1].index;
            }
        }
        delete[] sort_list;
    }

    void PLVF::UpdateNeighbours() {
        DistanceInfo *sort_list = new DistanceInfo[weight_num_];
        for (int i = 0; i < weight_num_; ++i) {
            for (int j = 0; j < weight_num_; ++j) {
                /**
                 * find neighbors by cosine, not euclidean distance
                 */
                sort_list[j].distance = CalculateCos(lambda_[i], lambda_[j], g_GlobalSettings->obj_num_);
                sort_list[j].index = j;

            }
            std::sort(sort_list, sort_list + weight_num_, [](const DistanceInfo &left, const DistanceInfo &right) {
                return left.distance > right.distance;
            });
            for (int j = 0; j < neighbour_num_; j++) {
                neighbour_[i][j] = sort_list[j + 1].index;
            }
        }
        delete[] sort_list;
    }

    void PLVF::Crossover(Individual **parent_pop, int current_index, Individual *offspring) {
        // randomly select two parent from current individual's neighbours
        int k = rnd(0, neighbour_num_ - 1);
        int l = rnd(0, neighbour_num_ - 1);
        Individual *parent1 = parent_pop[neighbour_[current_index][k]];
        Individual *parent2 = parent_pop[neighbour_[current_index][l]];

        SBX(parent1, parent2, g_GlobalSettings->offspring_population_[1], offspring, g_GlobalSettings);

    }

    void PLVF::UpdateSubproblem(Individual *offspring, int current_index, int aggregation_type) {
        double *offspring_fitness = new double[neighbour_num_];
        double *neighbour_fitness = new double[neighbour_num_];

        // calculate fitness;
        //default: aggregation_type := 0
        switch (aggregation_type) {
            case 0:
                // inverse chebycheff
                for (int i = 0; i < neighbour_num_; ++i) {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = CalInverseChebycheff(current_ind, lambda_[weight_index], ideal_point_,
                                                                g_GlobalSettings->obj_num_);
                    offspring_fitness[i] = CalInverseChebycheff(offspring, lambda_[weight_index], ideal_point_,
                                                                g_GlobalSettings->obj_num_);
                }
                break;

            case 1:
                // weighted sum
                for (int i = 0; i < neighbour_num_; ++i) {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = CalWeightedSum(current_ind, lambda_[weight_index], ideal_point_,
                                                          g_GlobalSettings->obj_num_);
                    offspring_fitness[i] = CalWeightedSum(offspring, lambda_[weight_index], ideal_point_,
                                                          g_GlobalSettings->obj_num_);
                }
                break;

            case 2:
                // PBI
                for (int i = 0; i < neighbour_num_; ++i) {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = CalPBI(current_ind, lambda_[weight_index], ideal_point_,
                                                  g_GlobalSettings->obj_num_, pbi_theta_);
                    offspring_fitness[i] = CalPBI(offspring, lambda_[weight_index], ideal_point_,
                                                  g_GlobalSettings->obj_num_, pbi_theta_);
                }
                break;
            case 3:
                // L_alpha_norm
                for (int i = 0; i < neighbour_num_; ++i) {
                    int weight_index = neighbour_[current_index][i];
                    Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
                    neighbour_fitness[i] = L_alpha_norm(lambda_[weight_index], current_ind->obj_, ideal_point_,
                                                        g_GlobalSettings->obj_num_, alpha_);
                    offspring_fitness[i] = L_alpha_norm(lambda_[weight_index], offspring->obj_, ideal_point_,
                                                        g_GlobalSettings->obj_num_, alpha_);
                }
                break;
            default:
                break;
        }

        // update subproblem
        for (int i = 0; i < neighbour_num_; ++i) {
            if (offspring_fitness[i] < neighbour_fitness[i]) {
                CopyIndividual(offspring, g_GlobalSettings->parent_population_[neighbour_[current_index][i]]);
            }
        }
        delete[] neighbour_fitness;
        delete[] offspring_fitness;
    }
}