#include <ctime>
#include <cstdio>
#include <pthread.h>
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <time.h>
#include <Python.h>

#include "core/file.h"
#include "core/global.h"
#include "core/individual.h"
#include "problem/zdt.h"
#include "problem/dtlz.h"
#include "algorithms/moead/moead.h"
#include "algorithms/moead_de/moead_de.h"
#include "algorithms/nemo_0/nemo_0.h"
#include "algorithms/iemod/iemod.h"
#include "algorithms/nsga2/nsga2.h" 
#include "algorithms/ibea/ibea.h"
#include "algorithms/spea2/spea2.h"
#include "algorithms/smsemoa/smsemoa.h"
#include "algorithms/moead_dra/moead_dra.h"
#include "algorithms/moead_frrmab/moead_frrmab.h"
#include "algorithms/hype/hype.h"
#include "metric/hv.h"
#include "metric/igd.h"
#include "metric/gd.h"
#include "metric/spacing.h"
#include "random/random.h"

#include "problem/correlated_500_knapsack.h"

#if defined(__linux) || defined(linux)
#include <sys/time.h>
#endif

using emoc::g_GlobalSettingsArray;
using emoc::EMOCParameters;


struct ThreadParamters
{
	EMOCParameters *para;

	int run_start;
	int run_end;
	int thread_id;
};


void *Work(void *args);
void EMOCMultiThreadTest(EMOCParameters *parameter);
void EMOCSingleThreadTest(EMOCParameters *parameter);

int main(int argc, char* argv[])
{
	// initilize some bases for random number
	randomize();
	// initialize parameters
	emoc::EMOCParameters *parameter = new emoc::EMOCParameters();
	ParseParamerters(argc, argv, parameter); // set params from file and command

	//ReadParametersFromFile("src/config/config.txt", parameter);
	//parameter->igd_value = (double *)malloc(sizeof(double) * parameter->runs_num);
    time_t now_time = time(NULL);
    tm* t_tm = localtime(&now_time);
    std::cout<<"local time: "<<asctime(t_tm)<<std::endl;
	std::cout << "current task:" << std::endl;
	std::cout << "-------------------------------------------" << std::endl;
	std::cout << "problem:              " << parameter->problem_name << std::endl;
	std::cout << "algorithm:            " << parameter->algorithm_name << std::endl;
	std::cout << "population number:    " << parameter->population_num << std::endl;
	std::cout << "decision number:      " << parameter->decision_num << std::endl;
	std::cout << "objective number:     " << parameter->objective_num << std::endl;
    std::cout << "DM's weight:          " << parameter->weight_stringType << std::endl;
	std::cout << "evaluation:           " << parameter->max_evaluation << std::endl;
    if (parameter->tau == 0)
    {
        std::cout << "ask time:             " << parameter->askTime << std::endl;
    }
    else
    {
        std::cout << "elicitation interval: " << parameter->tau << std::endl;
    }
    std::cout << "step size:            " << parameter->step_size << std::endl;
    std::cout << "inquiries number:     " << parameter->inquiriesNum << std::endl;
    std::cout << "runs:                 " << parameter->runs_num << std::endl;
	std::cout << "is open multithread:  " << parameter->is_open_multithread << std::endl;
	std::cout << "multithread number:   " << parameter->thread_num << std::endl;
    std::cout << "noise level:          " << parameter->kappa << std::endl;
	std::cout << "-------------------------------------------\n" << std::endl;

	clock_t start, end;
	start = clock();

    // read knapsack configuration from local file
    if (parameter->problem_name == "knapsack")
    {
        assert(parameter->objective_num<=10);
        using emoc::weight_knapsack;
        using emoc::profit_knapsack;
        using emoc::capacity_knapsack;
        char weight_file_name[256];
        char profit_file_name[256];
        sprintf(weight_file_name,"./configFiles/test_problem/knapsack/knapsack.500.%d.weight.txt", parameter->objective_num);
        sprintf(profit_file_name,"./configFiles/test_problem/knapsack/knapsack.500.%d.profit.txt", parameter->objective_num);
        FILE *fpt_weight = fopen(weight_file_name, "r");
        FILE *fpt_profit = fopen(profit_file_name, "r");
        std::ifstream in_weight, in_profit;
        in_weight.open(weight_file_name);
        in_profit.open(profit_file_name);
        int temp = 0;
        for (int i = 0; i < parameter->objective_num; i++)
        {
            temp = 0;
            for (int j = 0; j < 500; j++)
            {
                in_weight>>weight_knapsack[i][j];
                in_profit>>profit_knapsack[i][j];
                temp += weight_knapsack[i][j];
            }
            capacity_knapsack[i] = temp / 2;
        }
        in_weight.close();
        in_profit.close();
        #ifdef DEBUG
        std::cout<<"---------weight for knapsack---------"<<std::endl;
        for (int i = 0; i < parameter->objective_num; ++i)
        {
            for (int j = 0; j < 500; ++j)
            {
                std::cout<<weight_knapsack[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"---------profit for knapsack---------"<<std::endl;
        for (int i = 0; i < parameter->objective_num; ++i)
        {
            for (int j = 0; j < 500; ++j)
            {
                std::cout<<profit_knapsack[i][j]<<" ";
            }
            std::cout<<std::endl;
        }
        std::cout<<"---------capacity for knapsack---------"<<std::endl;
        for (int i = 0; i < parameter->objective_num; ++i)
        {
            std::cout<<capacity_knapsack[i]<<" ";
        }
        std::cout<<std::endl;
        #endif
    }

	// EMOC test run
    printf("<info> is_open_multithread: %d\n", parameter->is_open_multithread);
	if (parameter->is_open_multithread)
    {
        EMOCMultiThreadTest(parameter);
    }
	else
    {
        EMOCSingleThreadTest(parameter);
    }
	end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	printf("------total run time: %fs--------\n", time);

	/*for (int i = 0; i < parameter->runs_num; ++i)
	{
		printf("run %d igd value: %f \n", i, parameter->igd_value[i]);
	}*/

	delete parameter;
    std::cout<<"The run time is: "<<(double)clock() / CLOCKS_PER_SEC<<std::endl;
	/*std::cout<<"\nTask has finished, please enter to exit."<<std::endl;
	std::cin.get();*/
	
	return 0;
}

void *Work(void *args)
{
	ThreadParamters *parameter = (ThreadParamters *)args;
    const char *algorithm_name = parameter->para->algorithm_name.c_str();
    const char *problem_name = parameter->para->problem_name.c_str();
    int population_num = parameter->para->population_num;
    int dec_num = parameter->para->decision_num;
    int obj_num = parameter->para->objective_num;
    int max_eval = parameter->para->max_evaluation;
    int output_interval = parameter->para->output_interval;
    int askTime = parameter->para->askTime;
    double step_size = parameter->para->step_size;
    int tau = parameter->para->tau;
    int inquiriesNum= parameter->para->inquiriesNum;
    const char *weight_stringType = parameter->para->weight_stringType.c_str();
    double kappa=parameter->para->kappa;
	for (int run = parameter->run_start; run <= parameter->run_end; ++run)
	{
		int thread_id = parameter->thread_id;
        std::cout << "[INFO] Thread ID: " << thread_id << std::endl;
        std::cout << "[INFO] Run No: " << run << std::endl;
		// algorithm main entity
        g_GlobalSettingsArray[thread_id] = new emoc::Global(algorithm_name, problem_name, population_num, dec_num, obj_num, max_eval, thread_id, output_interval, run, askTime, tau, step_size, inquiriesNum, weight_stringType, kappa);
        g_GlobalSettingsArray[thread_id]->Start();
		/*std::string problem_name = g_GlobalSettingsArray[thread_id]->problem_name_;
		int obj_num = g_GlobalSettingsArray[thread_id]->obj_num_;
		double igd = emoc::CalculateIGD(g_GlobalSettingsArray[thread_id]->parent_population_.data(), g_GlobalSettingsArray[thread_id]->population_num_, obj_num, problem_name);
		printf("current thread id : %d, runs: %d, igd:%f\n", thread_id, run, igd);
		parameter->para->igd_value[run] = igd;*/
		//RecordPop(run, 0, g_GlobalSettingsArray[thread_id]);
		delete g_GlobalSettingsArray[thread_id];
	}
	return nullptr;
}

void EMOCMultiThreadTest(EMOCParameters *parameter)
{
	int thread_num = parameter->thread_num;
	std::vector<ThreadParamters*> thread_para_array(thread_num, nullptr);
	for (int i = 0; i < thread_num; ++i)
	{
		thread_para_array[i] = (ThreadParamters *)malloc(sizeof(ThreadParamters));
		thread_para_array[i]->para = parameter;
	}

	std::vector<int> job_overload(thread_num, 0);
	int interval = (double)parameter->runs_num / thread_num;
	int remainder = parameter->runs_num % thread_num;
	for (int i = 0; i < thread_num; ++i)
	{
		job_overload[i] = interval;
		if (remainder-- > 0)
        {
            job_overload[i]++;
        }
		//printf("thread %d: %d runs\n",i, job_overload[i]);
	}
	// multithread running
	std::vector<pthread_t> tid(thread_num);
	int total_overload = parameter->start_index;

    if (parameter->algorithm_name == "nsga2_ltr"
        || parameter->algorithm_name == "moead_ltr"
        || parameter->algorithm_name == "r2_ibea_ltr") {
        Py_Initialize(); // 初始化
        //判断初始化是否成功
        if(!Py_IsInitialized())
        {
            std::cerr<<"[ERROR] Python init failed!\n";
            return;
        }
        Py_DECREF(PyImport_ImportModule("threading"));//TODO: DEBUG
        PyEval_InitThreads(); // 初始化和获取全局解释器锁,初始化线程支持
        PyEval_ReleaseThread(PyThreadState_Get()); // 重置当前线程状态为NULL并释放全局解释器锁
    }

    // create child threads
	for (int i = 0; i < thread_num; ++i)
	{
		if (job_overload[i] > 0)
		{
			thread_para_array[i]->run_start = total_overload;
			thread_para_array[i]->run_end = total_overload + job_overload[i] - 1;
			thread_para_array[i]->thread_id = i;
			total_overload += job_overload[i];
		}
		else
        {
            continue;
        }
		pthread_create(&tid[i], nullptr, Work, (void *)thread_para_array[i]);
	}

	// multi-thread destroy
	for (int i = 0; i < thread_num; ++i)
	{
		if (job_overload[i] > 0)
        {
            pthread_join(tid[i], nullptr);
        }
	}

	for (int i = 0; i < thread_num; ++i)
    {
        free(thread_para_array[i]);
    }

    if (parameter->algorithm_name == "nsga2_ltr"
        || parameter->algorithm_name == "moead_ltr"
        || parameter->algorithm_name == "r2_ibea_ltr") {
        PyGILState_Ensure();
        Py_Finalize();
    }
}

void EMOCSingleThreadTest(EMOCParameters *parameter)
{
    Py_Initialize(); // 初始化
    //判断初始化是否成功
    if(!Py_IsInitialized())
    {
        std::cerr<<"[ERROR] Python init failed!\n";
        return;
    }

	const char *algorithm_name = parameter->algorithm_name.c_str();
	const char *problem_name = parameter->problem_name.c_str();
	int population_num = parameter->population_num;
	int dec_num = parameter->decision_num;
	int obj_num = parameter->objective_num;
	int max_eval = parameter->max_evaluation;
	int output_interval = parameter->output_interval;
	int askTime = parameter->askTime;
	double step_size = parameter->step_size;
	int tau = parameter->tau;
	int inquiriesNum= parameter->inquiriesNum;
    const char *weight_stringType = parameter->weight_stringType.c_str();
    double kappa=parameter->kappa;


    for (int run = parameter->start_index; run < parameter->start_index + parameter->runs_num; ++run)
	{
		int thread_id = 0;
        std::cout << "[INFO] Thread ID: " << thread_id << std::endl;
        std::cout << "[INFO] Run No: " << run << std::endl;
		//run time recording
		clock_t start, end;
		start = clock();

		// algorithm main entity
        g_GlobalSettingsArray[thread_id] = new emoc::Global(algorithm_name, problem_name, population_num, dec_num, obj_num, max_eval, thread_id, output_interval, run, askTime, tau, step_size, inquiriesNum, weight_stringType, kappa);
		g_GlobalSettingsArray[thread_id]->Start();
        end = clock();
		double time = (double)(end - start) / CLOCKS_PER_SEC;

		printf("run %d time: %fs\n",run, time);
		delete g_GlobalSettingsArray[thread_id];
	}
    Py_Finalize();
}