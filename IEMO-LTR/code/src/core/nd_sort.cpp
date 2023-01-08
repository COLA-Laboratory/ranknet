#include "core/nd_sort.h"

#include <cstring>
#include <cstdlib>
#include <iostream>

#include "core/utility.h"

namespace emoc {

    void NonDominatedSort(Individual **pop, int pop_num, int obj_num)
    {
        int index = 0;
        int dominate_relation = 0;
        int current_rank = 0, unrank_num = pop_num;

        int *ni = nullptr;             // store the number of points that dominate i-th solution
        int **si = nullptr;            // store the solution index of which i-th solution dominates
        int	*Q = nullptr;              // store the solution which ni is 0
        int *dominate_num = nullptr;   // store the number of dominate points of i-th solution
        Individual *ind_tempA = nullptr, *ind_tempB = nullptr;

        ni = (int *)malloc(sizeof(int) * pop_num);
        memset(ni, 0, sizeof(int) * pop_num);

        si = (int **)malloc(sizeof(int *) * pop_num);
        for (int i = 0; i < pop_num; i++)
        {
            si[i] = (int *)malloc(sizeof(int) * pop_num);
            memset(si[i], 0, sizeof(int) * pop_num);
        }

        Q = (int *)malloc(sizeof(int) * pop_num);
        memset(Q, 0, sizeof(int) * pop_num);

        dominate_num = (int*)malloc(sizeof(int) * pop_num);
        memset(dominate_num, 0, sizeof(int) * pop_num);

        // FIXME
        int domination_matrix[pop_num][pop_num];
        memset(domination_matrix, 0, sizeof(int) * pop_num * pop_num);
        for (int i = 0; i < pop_num - 1; i++)
        {
            for (int j = i + 1; j < pop_num; ++j)
            {
                domination_matrix[i][j] = CheckDominance(pop[i], pop[j], obj_num);
                domination_matrix[j][i] = -domination_matrix[i][j];
            }
        }

        for (int i = 0; i < pop_num; ++i)
        {
            index = 0;
            for (int j = 0; j < pop_num; ++j)
            {
                if (DOMINATE == domination_matrix[i][j])
                {
                    si[i][index++] = j;
                }
                else if (DOMINATED == domination_matrix[i][j])
                {
                    ++ni[i];
                }
            }
            // number of solutions that solution i dominates
            dominate_num[i] = index;
        }


//		for (int i = 0; i < pop_num; i++)
//		{
//			ind_tempA = pop[i];
//			index = 0;
//			for (int j = 0; j < pop_num; j++)
//			{
//				if (i == j)
//                {
//                    continue;
//                }
//
//				ind_tempB = pop[j];
//
//				dominate_relation = check_r_dominance(ind_tempA, ind_tempB, obj_num);
//                if (DOMINATE == dominate_relation)
//				{
//					si[i][index++] = j;
//				}
//				else if (DOMINATED == dominate_relation)
//				{
//					++ni[i];
//				}
//				else;
//			}
//            // number of solutions that solution i dominates
//			dominate_num[i] = index;
//		}



        // assign domination rank
        bool flag_test = false;
        while (unrank_num)
        {
            // std::cerr<<current_rank<<std::endl;
            index = 0;
            flag_test = false;
            for (int i = 0; i < pop_num; i++)
            {
                if (ni[i] == 0)
                {// no solution dominates solution i
                    flag_test = true;
                    pop[i]->rank_ = current_rank;
                    Q[index++] = i;
                    --unrank_num;
                    ni[i] = -1;
                }
            }
            if (!flag_test)
            {
                for (int ii = 0; ii < pop_num; ii++)
                {
                    std::cerr<<ni[ii]<<" ";
                }
                std::cerr<<std::endl;
                exit(-1);
            }

            ++current_rank;

            for (int i = 0; i < index; i++)
            {
                for (int j = 0; j < dominate_num[Q[i]]; j++)
                {
                    // Q[i] 表示支配它的解的数量为0的解的下标
                    // si[Q[i]][j] 表示下标为Q[i]的解支配的解的下标
                    // ni表示支配这个解的解个数
                    --ni[si[Q[i]][j]];
                }
            }
        }

        free(ni);
        for (int i = 0; i < pop_num; i++)
        {
            free(si[i]);
        }
        free(si);
        free(Q);
        free(dominate_num);
        return;
    }

	void rNonDominatedSort(Individual **pop, int pop_num, int obj_num)
	{
		int index = 0; 
		int dominate_relation = 0;
		int current_rank = 0, unrank_num = pop_num; 

		int *ni = nullptr;             // store the number of points that dominate i-th solution
		int **si = nullptr;            // store the solution index of which i-th solution dominates
		int	*Q = nullptr;              // store the solution which ni is 0
		int *dominate_num = nullptr;   // store the number of dominate points of i-th solution
		Individual *ind_tempA = nullptr, *ind_tempB = nullptr;

		ni = (int *)malloc(sizeof(int) * pop_num);
		memset(ni, 0, sizeof(int) * pop_num);

		si = (int **)malloc(sizeof(int *) * pop_num);
		for (int i = 0; i < pop_num; i++)
		{
			si[i] = (int *)malloc(sizeof(int) * pop_num);
			memset(si[i], 0, sizeof(int) * pop_num);
		}

		Q = (int *)malloc(sizeof(int) * pop_num);
		memset(Q, 0, sizeof(int) * pop_num);

		dominate_num = (int*)malloc(sizeof(int) * pop_num);
		memset(dominate_num, 0, sizeof(int) * pop_num);

        // FIXME
        int domination_matrix[pop_num][pop_num];
        memset(domination_matrix, 0, sizeof(int) * pop_num * pop_num);
        for (int i = 0; i < pop_num - 1; i++)
        {
            for (int j = i + 1; j < pop_num; ++j)
            {
                domination_matrix[i][j] = check_r_dominance(pop[i], pop[j], obj_num);
                domination_matrix[j][i] = -domination_matrix[i][j];
            }
        }

        for (int i = 0; i < pop_num; ++i)
        {
            index = 0;
            for (int j = 0; j < pop_num; ++j)
            {
                if (DOMINATE == domination_matrix[i][j])
                {
                    si[i][index++] = j;
                }
                else if (DOMINATED == domination_matrix[i][j])
                {
                    ++ni[i];
                }
            }
            // number of solutions that solution i dominates
            dominate_num[i] = index;
        }


//		for (int i = 0; i < pop_num; i++)
//		{
//			ind_tempA = pop[i];
//			index = 0;
//			for (int j = 0; j < pop_num; j++)
//			{
//				if (i == j)
//                {
//                    continue;
//                }
//
//				ind_tempB = pop[j];
//
//				dominate_relation = check_r_dominance(ind_tempA, ind_tempB, obj_num);
//                if (DOMINATE == dominate_relation)
//				{
//					si[i][index++] = j;
//				}
//				else if (DOMINATED == dominate_relation)
//				{
//					++ni[i];
//				}
//				else;
//			}
//            // number of solutions that solution i dominates
//			dominate_num[i] = index;
//		}



        // assign domination rank
        bool flag_test = false;
		while (unrank_num)
		{
            // std::cerr<<current_rank<<std::endl;
			index = 0;
            flag_test = false;
			for (int i = 0; i < pop_num; i++)
			{
				if (ni[i] == 0)
				{// no solution dominates solution i
                    flag_test = true;
					pop[i]->rank_ = current_rank;
					Q[index++] = i;
					--unrank_num;
					ni[i] = -1;
				}
			}
            if (!flag_test)
            {
                for (int ii = 0; ii < pop_num; ii++)
                {
                    std::cerr<<ni[ii]<<" ";
                }
                std::cerr<<std::endl;
                exit(-1);
            }

			++current_rank;

			for (int i = 0; i < index; i++)
			{
				for (int j = 0; j < dominate_num[Q[i]]; j++)
				{
                    // Q[i] 表示支配它的解的数量为0的解的下标
                    // si[Q[i]][j] 表示下标为Q[i]的解支配的解的下标
                    // ni表示支配这个解的解个数
					--ni[si[Q[i]][j]];
				}
			}
		}

		free(ni);
		for (int i = 0; i < pop_num; i++)
		{
			free(si[i]);
		}
		free(si);
		free(Q);
		free(dominate_num);
        std::cerr<<"rNonDominatedSort() returns."<<std::endl;
        return;
	}

}