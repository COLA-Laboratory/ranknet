#include "algorithms/moead_dra/moead_dra.h"

#include <cmath>
#include <algorithm>

#include "core/global.h"
#include "core/utility.h"
#include "core/uniform_point.h"
#include "operator/de.h"
#include "operator/mutation.h"
#include "random/random.h"

namespace emoc {

	MOEADDRA::MOEADDRA(Problem *problem, int thread_num) :
		Algorithm(problem, thread_num),
		lambda_(nullptr),
		weight_num_(0),
		neighbour_(nullptr),
		neighbour_selectpro_(0.9),
		ideal_point_(new double[g_GlobalSettings->obj_num_])
	{

	}

	MOEADDRA::~MOEADDRA()
	{
		for (int i = 0; i < weight_num_; ++i)
		{
			delete[] lambda_[i];
			delete[] neighbour_[i];
			lambda_[i] = nullptr;
			neighbour_[i] = nullptr;
		}
		delete[] lambda_;
		delete[] neighbour_;
		delete[] ideal_point_;
		delete[] selected_indices_;
		delete[] old_obj_;
		delete[] utility_;
		delete[] delta_;
		lambda_ = nullptr;
		neighbour_ = nullptr;
		ideal_point_ = nullptr;
		selected_indices_ = nullptr;
		old_obj_ = nullptr;
		utility_ = nullptr;
		delta_ = nullptr;
	}

	void MOEADDRA::Run()
	{
		Initialization();

		// calculate and store the old fitness
		CalculateFitness(g_GlobalSettings->parent_population_.data(), weight_num_, old_obj_);

		Individual *offspring = g_GlobalSettings->offspring_population_[0];
		while (!g_GlobalSettings->IsTermination())
		{
			// begin each iteration
			g_GlobalSettings->iteration_num_++;
			for (int subgeneration = 0; subgeneration < 5; ++subgeneration)
			{
				SelectCurrentSubproblem();
				for (int i = 0; i < selected_size_; ++i)
				{
					int index = selected_indices_[i];
					// set current iteration's neighbour type
					if (randomperc() < neighbour_selectpro_)
						neighbour_type_ = NEIGHBOUR;
					else
						neighbour_type_ = GLOBAL;

					// generate offspring for current subproblem
					Crossover(g_GlobalSettings->parent_population_.data(), index, offspring);
					MutationInd(offspring, g_GlobalSettings);
					EvaluateInd(offspring);

					// update ideal point
					UpdateIdealpoint(offspring, ideal_point_, g_GlobalSettings->obj_num_);

					// update neighbours' subproblem 
					UpdateSubproblem(offspring, index);
				}
			}

			if (g_GlobalSettings->iteration_num_ % 10 == 0)
			{
				UpdateUtility();
			}	

			// record the population every interval generations and the first and last genration 
			if (g_GlobalSettings->iteration_num_ % g_GlobalSettings->output_interval_ == 0 || g_GlobalSettings->iteration_num_ == 1
				|| g_GlobalSettings->IsTermination())
				TrackPopulation(g_GlobalSettings->iteration_num_);
		}
	}

	void MOEADDRA::Initialization()
	{
		// initialize parent population
		g_GlobalSettings->InitializePopulation(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);
		EvaluatePop(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);

		// generate weight vectors
		lambda_ = UniformPoint(g_GlobalSettings->population_num_, &weight_num_, g_GlobalSettings->obj_num_);
		replace_num_ = (weight_num_ / 100) ? (weight_num_ / 100) : 2;

		// set the neighbours of each individual
		SetNeighbours();

		// initialize ideal point
		UpdateIdealpoint(g_GlobalSettings->parent_population_.data(), weight_num_, ideal_point_, g_GlobalSettings->obj_num_);

		// set selected size
		selected_size_ = weight_num_ / 5;
		selected_indices_ = new int[selected_size_];

		// initialize utility related data
		old_obj_ = new double[weight_num_];
		delta_ = new double[weight_num_];
		utility_ = new double[weight_num_];
		for (int i = 0; i < weight_num_; ++i)
			utility_[i] = 1.0;
	}

	void MOEADDRA::SetNeighbours()
	{
		// set neighbour size and allocate memory
		neighbour_num_ = weight_num_ / 10;
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

			std::sort(sort_list, sort_list + weight_num_, [](const DistanceInfo &left, const DistanceInfo &right) {
				return left.distance < right.distance;
			});

			for (int j = 0; j < neighbour_num_; j++)
			{
				neighbour_[i][j] = sort_list[j + 1].index;
			}
		}

		delete[] sort_list;
	}

	void MOEADDRA::Crossover(Individual **parent_pop, int current_index, Individual *offspring)
	{
		int size = neighbour_type_ == NEIGHBOUR ? neighbour_num_ : weight_num_;
		int parent2_index = 0, parent3_index = 0;

		// randomly select two parents according to neighbour type
		if (neighbour_type_ == NEIGHBOUR)
		{
			parent2_index = neighbour_[current_index][rnd(0, size - 1)];
			parent3_index = neighbour_[current_index][rnd(0, size - 1)];
		}
		else
		{
			parent2_index = rnd(0, size - 1);
			parent3_index = rnd(0, size - 1);
		}

		Individual *parent1 = parent_pop[current_index];
		Individual *parent2 = parent_pop[parent2_index];
		Individual *parent3 = parent_pop[parent3_index];
		DE(parent1, parent2, parent3, offspring, g_GlobalSettings);
	}

	void MOEADDRA::UpdateSubproblem(Individual *offspring, int current_index)
	{
		int size = neighbour_type_ == NEIGHBOUR ? neighbour_num_ : weight_num_;
		int *perm_index = new int[size];
		random_permutation(perm_index, size);

		int count = 0, weight_index = 0;
		double offspring_fitness = 0.0;
		double neighbour_fitness = 0.0;

		// calculate fitness and update subproblem;
		for (int i = 0; i < size; ++i)
		{
			if (count >= replace_num_)
				break;

			if (neighbour_type_ == NEIGHBOUR)
				weight_index = neighbour_[current_index][perm_index[i]];
			else
				weight_index = perm_index[i];

			Individual *current_ind = g_GlobalSettings->parent_population_[weight_index];
			offspring_fitness = CalInverseChebycheff(offspring, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_);
			neighbour_fitness = CalInverseChebycheff(current_ind, lambda_[weight_index], ideal_point_, g_GlobalSettings->obj_num_);
			if (offspring_fitness < neighbour_fitness)
			{
				CopyIndividual(offspring, g_GlobalSettings->parent_population_[weight_index]);
				count++;
			}
		}

		delete[] perm_index;
	}

	void MOEADDRA::SelectCurrentSubproblem()
	{
		int rand[10] = { 0 }, current_max_index = 0;

		// add boundry first
		for (int i = 0; i < g_GlobalSettings->obj_num_; ++i)
		{
			for (int j = 0; j < weight_num_; ++j)
			{
				if (fabs(lambda_[j][i] - 1) < EPS)
					selected_indices_[i] = j;
			}
		}

		// 10-tournment selection
		for (int i = g_GlobalSettings->obj_num_; i < selected_size_; ++i)
		{
			for (int j = 0; j < 10; ++j)
			{
				rand[j] = rnd(0, weight_num_ - 1);
			}

			double max = -INF;
			for (int j = 0; j < 10; ++j)
			{
				if (utility_[rand[j]] > max)
				{
					max = utility_[rand[j]];
					current_max_index = rand[j];
				}
			}
			selected_indices_[i] = current_max_index;
		}
	}

	void MOEADDRA::CalculateFitness(Individual **pop, int pop_num, double *fitness)
	{
		for (int i = 0; i < pop_num; ++i)
		{
			double fit = CalInverseChebycheff(pop[i], lambda_[i], ideal_point_, g_GlobalSettings->obj_num_);
			fitness[i] = fit;
		}
	}

	void MOEADDRA::UpdateUtility()
	{
		// calculate delta
		for (int i = 0; i < weight_num_; ++i)
		{
			delta_[i] = fabs(g_GlobalSettings->parent_population_[i]->fitness_ - old_obj_[i]) / old_obj_[i];
			old_obj_[i] = g_GlobalSettings->parent_population_[i]->fitness_;
		}

		// update utility
		for (int i = 0; i < weight_num_; ++i)
		{
			if (delta_[i] > 0.001)
			{
				utility_[i] = 1.0;
			}
			else
			{
				utility_[i] = utility_[i] * (0.95 + 0.05 * (delta_[i] / 0.001));
			}
		}
	}
}