// This file provides some important utility functions
#pragma once
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <numeric>

#include "core/global.h"
#include "core/individual.h"
#include <iostream>
#include "eigen3/Eigen/SVD"
namespace emoc {

	enum DominateReleation 
	{
		DOMINATED = -1,
		NON_DOMINATED = 0,
		DOMINATE = 1
	};
	DominateReleation CheckDominance(Individual *ind1, Individual *ind2, int obj_num);
    DominateReleation check_r_dominance(Individual *ind1, Individual *ind2, int obj_num);
	int WeaklyDominates(double *point1, double *point2, int obj_num);

	int Combination(int n, int k);
	double CalculateDotProduct(double *vector1, double *vector2, int dimension);
	double CalculateCos(double *a, double *b, int dimension);
    double CalculateCos(double *a, double *b, int dimension, int power);
    double L_alpha_norm(double *weight, double *obj, double *referencePoint, int dimension, double alpha);
    double L_alpha_norm_DM(double *weight, double *obj, int dimension, double alpha);
	double CalculateSin(double *a, double *b, int dimension);
	double CalculateNorm(double *vector, int dimension);
	double CalEuclidianDistance(double *a, double *b, int dimension);
	double CalPerpendicularDistance(double *a, double *weight, int dimension);

	void UpdateIdealpoint(Individual *ind, double *ideal_point, int obj_num);
	void UpdateNadirpoint(Individual *ind, double *nadir_point, int obj_num);
	void UpdateIdealpoint(Individual **pop, int pop_num, double *ideal_point, int obj_num);
	void UpdateNadirpoint(Individual **pop, int pop_num, double *nadir_point, int obj_num);

	// aggregation functions
	double CalWeightedSum(Individual *ind, double *weight_vector, double *ideal_point, int obj_num);
	double CalInverseChebycheff(Individual *ind, double *weight_vector, double *ideal_point, int obj_num);
    double CalInverseChebycheff(double *obj, double *weight_vector, double *ideal_point, int obj_num);
    double CalInverseChebycheff(Individual *ind, double *weight_vector, int obj_num);
    double CalInverseChebycheff(double *obj, double *weight_vector, int obj_num);
    double CalChebycheff(Individual *ind, double *weight_vector, int obj_num);
	double CalPBI(Individual *ind, double *weight_vector, double *ideal_point, int obj_num, double theta = 0.0 );

    double** calculatePinv(double** n, int row, int col);
    void TrainRBFNet(double **c, double *weight,double *out,  double sigma, int size, int dimension);
    void UsingRBFNet(Individual *ind, double **c, double *weight, double sigma, int size, int dimension);

    double NDCG();
    double topK(double *model, double *gold, int len, int k);
    double* SetWeight(const std::string& weightstring);

    void InitializeReferencePoint(double *reference_point, int obj_num);
    void UpdateReferencePoint(Individual **mixed_pop, int mixed_pop_num, int obj_num, double *reference_point);
    void UpdateIdealPoint(Individual **pop, int pop_num, int obj_num, double *ideal_point);
    void display_pop(FILE *gp,Individual **pop,int pop_num, int obj_num, int gen);
    void display_mixed_pop(FILE *gp, int obj_num, int gen);
    void GNUPlot_score(FILE *gp,Individual **pop,int pop_num, int obj_num, int gen);
    /**
     * convert real value in [0,1] to 0 or 1
     * @param pop
     * @param pop_num
     */
    void initialize_bi_population(Individual **pop, int pop_num, int dec_num, int obj_num);
    extern int weight_knapsack[10][500];
    extern int profit_knapsack[10][500];
    extern int capacity_knapsack[10];
    bool isFeasible_knapsack(Individual *ind, int dec_num, int obj_num);
    void greedy_repair(Individual **pop, int pop_num, int dec_num, int obj_num);
    void greedy_repair_ind(Individual *ind, int dec_num, int obj_num);
    /**
     * get an integer in [min, max]
     * @param min
     * @param max
     * @return
     */
    int get_random_integer(int min, int max);


}