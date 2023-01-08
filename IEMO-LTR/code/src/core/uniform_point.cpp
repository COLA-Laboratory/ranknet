#include "core/uniform_point.h"
#include <cstdlib>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <pthread.h>
namespace emoc {

	static void SetWeight(double *weight, double unit, double sum, int dim, int *column, double **lambda, int obj_num)
	{
		int i;

		if (dim == obj_num)
		{
			for (i = 0; i < obj_num; i++)
				weight[i] = 0;
		}

		if (dim == 1)
		{
			weight[0] = unit - sum;
			for (i = 0; i < obj_num; i++)
				lambda[*column][i] = weight[i];
			*column = *column + 1;
			return;
		}
		for (i = 0; i <= unit - sum; i++)
		{
			weight[dim - 1] = i;
			SetWeight(weight, unit, sum + i, dim - 1, column, lambda, obj_num);
		}

		return;
	}

	double** UniformPoint(int num, int *weight_num, int obj_num)
	{
		
		int column = 0;
		double *vec = nullptr;
		double **lambda = nullptr;

		int gaps = 1;
		*weight_num = 0;
		while (1)
		{
			int layer_size = Combination(obj_num + gaps - 1, gaps);

			if (layer_size > num) break;
			*weight_num = layer_size;
			gaps = gaps + 1;
		}

		gaps = gaps - 1;
		lambda = new double*[*weight_num];
		for (int i = 0; i < *weight_num; i++)
		{
			lambda[i] = new double[obj_num];
		}

		vec = new double[obj_num];
		for (int i = 0; i < obj_num; i++)
			vec[i] = 0;
		SetWeight(vec, gaps, 0, obj_num, &column, lambda,obj_num);

		for (int i = 0; i < *weight_num; i++)
			for (int j = 0; j < obj_num; j++) 
				lambda[i][j] = lambda[i][j] / gaps;
			
		delete[] vec;
		return lambda;
	}

    double** HitAndRun(int num, int obj_num, char *run, char *outputFile)
    {
	    // signal(SIGCHLD, SIG_IGN);
        // usleep(100000);// sleep for 0.1 second
        pid_t status = system(run);
        /*int pid;
        #if DEBUG
        if (-1 == (pid=wait(NULL)))
        {
            std::cout<<"[INFO] No child process"<<std::endl;
        }
        else
        {
            std::cout<<"[INFO] Child process PID: "<<pid<<std::endl;
        }
        #endif*/
        if (-1 == status)
        {
            std::cerr<<"system error!\n";
            pthread_exit(0);
        }
        /*else
        {
            printf("exit status value = [0x%x]\n", status);
            if (WIFEXITED(status))
            {
                if (0 == WEXITSTATUS(status))
                {
                    printf("run shell script successfully.\n");
                }
                else
                {
                    printf("run shell script fail, script exit code: %d\n", WEXITSTATUS(status));
                }
            }
            else
            {
                printf("exit status = [%d]\n", WEXITSTATUS(status));
            }
        }*/
        double **lambda = nullptr;
        lambda = new double*[num];
        for (int i = 0; i < num; i++)
        {
            lambda[i] = new double[obj_num];
        }
        std::ifstream in;
        in.open(outputFile);
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < obj_num; j++)
            {
                in>>lambda[i][j];
            }
        }
        in.close();
        
        return lambda;
    }

    double** LoadUniformWeights(int num, int obj_num, char* file)
    {
        double **lambda = nullptr;
        lambda = new double*[num];
        for (int i = 0; i < num; i++)
        {
            lambda[i] = new double[obj_num];
        }
        std::ifstream in;
        in.open(file);
        for (int i = 0; i < num; i++)
        {
            for (int j = 0; j < obj_num; j++)
            {
                in>>lambda[i][j];
            }
        }
        in.close();
        return lambda;
    }



    double* GoldenPoint(std::string testProblem, int obj_num)
    {
        double *goldenPoint = nullptr;
        goldenPoint = new double[obj_num];
        if (testProblem == "dtlz1")
        {
            for (int i = 0; i < obj_num; ++i)
            {
                goldenPoint[i] = 0.5 / obj_num;
            }
            /*switch (obj_num)
            {
                case 3:
                    goldenPoint[0] = 0.2;
                    goldenPoint[1] = 0.15;
                    goldenPoint[2] = 0.15;
                    break;
                case 5:
                    goldenPoint[0] = 0.05;
                    goldenPoint[1] = 0.045;
                    goldenPoint[2] = 0.04;
                    goldenPoint[3] = 0.04;
                    goldenPoint[4] = 0.325;
                    break;
                case 8:
                    goldenPoint[0] = 0.25;
                    goldenPoint[1] = 0.035;
                    goldenPoint[2] = 0.040;
                    goldenPoint[3] = 0.035;
                    goldenPoint[4] = 0.035;
                    goldenPoint[5] = 0.030;
                    goldenPoint[6] = 0.040;
                    goldenPoint[7] = 0.035;
                    break;
                case 10:
                    goldenPoint[0] = 0.025;
                    goldenPoint[1] = 0.030;
                    goldenPoint[2] = 0.025;
                    goldenPoint[3] = 0.030;
                    goldenPoint[4] = 0.030;
                    goldenPoint[5] = 0.250;
                    goldenPoint[6] = 0.025;
                    goldenPoint[7] = 0.030;
                    goldenPoint[8] = 0.025;
                    goldenPoint[9] = 0.030;
                    break;
                default:
                    break;
            }*/

        }
        else if (testProblem == "dtlz2" || testProblem == "dtlz3" || testProblem == "dtlz4" || testProblem == "dtlz5" || testProblem == "dtlz6")
        {
            for (int i = 0; i < obj_num; ++i)
            {
                goldenPoint[i] = sqrt(1.0 / obj_num);
            }
            /*switch (obj_num)
            {
                case 3:
                    goldenPoint[0] = 0.686;
                    goldenPoint[1] = 0.514;
                    goldenPoint[2] = 0.514;
                    break;
                case 5:
                    goldenPoint[0] = 0.119;
                    goldenPoint[1] = 0.965;
                    goldenPoint[2] = 0.149;
                    goldenPoint[3] = 0.134;
                    goldenPoint[4] = 0.119;
                    break;
                case 8:
                    goldenPoint[0] = 0.131;
                    goldenPoint[1] = 0.112;
                    goldenPoint[2] = 0.150;
                    goldenPoint[3] = 0.131;
                    goldenPoint[4] = 0.935;
                    goldenPoint[5] = 0.131;
                    goldenPoint[6] = 0.150;
                    goldenPoint[7] = 0.131;
                    break;
                case 10:
                    goldenPoint[0] = 0.114;
                    goldenPoint[1] = 0.114;
                    goldenPoint[2] = 0.948;
                    goldenPoint[3] = 0.095;
                    goldenPoint[4] = 0.114;
                    goldenPoint[5] = 0.095;
                    goldenPoint[6] = 0.114;
                    goldenPoint[7] = 0.095;
                    goldenPoint[8] = 0.114;
                    goldenPoint[9] = 0.095;
                    break;
                default:
                    break;
            }*/
        }
        else if (testProblem == "dtlz3" || testProblem == "dtlz4")
        {
            /*switch (obj_num)
            {
                case 3:
                    goldenPoint[0] = 0.686;
                    goldenPoint[1] = 0.514;
                    goldenPoint[2] = 0.514;
                    break;
                case 5:
                    goldenPoint[0] = 0.445;
                    goldenPoint[1] = 0.4;
                    goldenPoint[2] = 0.533;
                    goldenPoint[3] = 0.4;
                    goldenPoint[4] = 0.445;
                    break;
                case 8:
                    goldenPoint[0] = 0.37;
                    goldenPoint[1] = 0.37;
                    goldenPoint[2] = 0.342;
                    goldenPoint[3] = 0.37;
                    goldenPoint[4] = 0.399;
                    goldenPoint[5] = 0.342;
                    goldenPoint[6] = 0.313;
                    goldenPoint[7] = 0.313;
                    break;
                case 10:
                    goldenPoint[0] = 0.345;
                    goldenPoint[1] = 0.377;
                    goldenPoint[2] = 0.345;
                    goldenPoint[3] = 0.251;
                    goldenPoint[4] = 0.314;
                    goldenPoint[5] = 0.283;
                    goldenPoint[6] = 0.314;
                    goldenPoint[7] = 0.283;
                    goldenPoint[8] = 0.345;
                    goldenPoint[9] = 0.283;
                    break;
                default:
                    break;
            }*/
        }
        else if (testProblem == "zdt1")
        {
            goldenPoint[0] = 0.3054;
            goldenPoint[1] = 0.4474;
        }
        return goldenPoint;
    }
}