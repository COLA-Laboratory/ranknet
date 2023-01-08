#define _CRT_SECURE_NO_WARNINGS

#include "core/file.h"
#include "utility.h"

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <fstream>

#if defined(_WIN32)
#include <direct.h>
#elif defined(__linux) || defined(linux)
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#endif

namespace emoc {

	static void SetParameter(char *parameter_name, char *value, EMOCParameters *para)
	{
		if (*parameter_name != '-')//first character
		{
			std::cout << "The parameter name should begin with '-'." << std::endl;
			std::cout << "Press enter to exit" << std::endl;
			std::cin.get();
			exit(-1);
		}

		parameter_name++;//skip '-'
		if (!strcmp(parameter_name, "algorithm"))
		{
			para->algorithm_name = value;
		}
		else if (!strcmp(parameter_name, "problem"))
		{
			para->problem_name = value;
		}
		else if (!strcmp(parameter_name, "D"))
		{
			para->decision_num = atoi(value);
		}
		else if (!strcmp(parameter_name, "M"))
		{
			para->objective_num = atoi(value);
		}
		else if (!strcmp(parameter_name, "N"))
		{
			para->population_num = atoi(value);
		}
		else if (!strcmp(parameter_name, "evaluation"))
		{
			para->max_evaluation = atoi(value);
		}
		else if (!strcmp(parameter_name, "save"))
		{
			para->output_interval = atoi(value);
		}
		else if (!strcmp(parameter_name, "run"))
		{
			para->runs_num = atoi(value);
		}
		else if (!strcmp(parameter_name, "thread"))
		{
            // 开启多线程
			para->is_open_multithread = 1;
			para->thread_num = atoi(value);
		}
        else if (!strcmp(parameter_name, "askTime"))
        {
            para->askTime = atoi(value);
        }
        else if (!strcmp(parameter_name, "stepSize"))
        {
            para->step_size = atof(value);
        }
        else if (!strcmp(parameter_name, "elicitationInterval"))
        {
            para->tau = atoi(value);
        }
        else if (!strcmp(parameter_name, "inquiriesNum"))
        {
            para->inquiriesNum = atoi(value);
        }
        else if (!strcmp(parameter_name, "W"))
        {
            para->weight_stringType = value;
        }
        else if (!strcmp(parameter_name, "start_index"))
        {
            para->start_index = atoi(value);
        }
        else if (!strcmp(parameter_name, "noise"))
        {
            para->kappa = atof(value);
        }
		else
		{
			std::cout << "Set a nonexistent parameter name." << std::endl;
			std::cout << "The parameters should be set like '-parameter_name value', e.g. '-algorithm nsga2'." << std::endl;
			std::cout << "Press enter to exit" << std::endl;
			std::cin.get();
			exit(-1);
		}
	}

	void PrintObjective(const char *filename, int dec_num, int obj_num, Individual** pop_table, int pop_num, const std::string &problem_name)
	{
        #if DEBUG
	    std::cout<<"[INFO] Output path/file: "<<filename<<std::endl;
        #endif

		FILE *fpt = fopen(filename, "w");
		if(fpt == nullptr)
		{
			std::cerr <<filename<< " doesn't exist." << std::endl;
			std::cerr << "Press enter to exit" << std::endl;
			std::cin.get();
			exit(-1);
		}
		for (int i = 0; i < pop_num; ++i)
		{
			Individual *ind = pop_table[i];

            // record the decision variables
            for (int j = 0; j < dec_num; ++j)
            {
                fprintf(fpt, "%lf,", ind->dec_[j]);
            }

            if (problem_name == "knapsack")
            {
                for (int j = 0; j < obj_num-1; ++j)
                {
                    fprintf(fpt, "%lf,", -1 * ind->obj_[j]);
                }
                fprintf(fpt, "%lf\n", -1 * ind->obj_[obj_num-1]);
            }
            else if(problem_name == "swimmer") {
                for (int j = 0; j < obj_num-1; ++j)
                {
                    fprintf(fpt, "%lf,", 300 - ind->obj_[j]);
                }
                fprintf(fpt, "%lf\n", 300 - ind->obj_[obj_num-1]);
            } else
            {
                for (int j = 0; j < obj_num-1; ++j)
                {
                    fprintf(fpt, "%lf,", ind->obj_[j]);
                }
                fprintf(fpt, "%lf\n", ind->obj_[obj_num-1]);
            }
		}
        if (fflush(fpt) == EOF )
        {
            std::cerr <<filename<< " flush error." << std::endl;
            std::cerr << "Press enter to exit" << std::endl;
            std::cin.get();
            exit(-1);
        }
		fclose(fpt);
	}
    void PrintObjective(const char *filename_pop, const char *filename_extremeSolution, const char *filename_extremeEP, int obj_num, Individual **pop_table, int pop_num, double *reference_point)
    {
        FILE *fpt1 = fopen(filename_pop, "w");
        FILE *fpt2 = fopen(filename_extremeSolution, "w");
        FILE *fpt3 = fopen(filename_extremeEP, "w");
        if(fpt1 == nullptr)
        {
            std::cout <<filename_pop<< " doesn't exist." << std::endl;
            std::cout << "Press enter to exit" << std::endl;
            std::cin.get();
            exit(-1);
        }
        if(fpt2 == nullptr)
        {
            std::cout <<filename_extremeSolution<< " doesn't exist." << std::endl;
            std::cout << "Press enter to exit" << std::endl;
            std::cin.get();
            exit(-1);
        }
        if(fpt3 == nullptr)
        {
            std::cout <<filename_extremeEP<< " doesn't exist." << std::endl;
            std::cout << "Press enter to exit" << std::endl;
            std::cin.get();
            exit(-1);
        }
        int bestIdx, worstIdx;
        double EP, bestEP, worstEP;
        bestEP = INF;
        worstEP = EPS;
        for (int i = 0; i < pop_num; ++i)
        {
            Individual *ind = pop_table[i];
            for (int j = 0; j < obj_num - 1; ++j)
            {
                fprintf(fpt1, "%f\t", ind->obj_[j]);
            }
            fprintf(fpt1, "%f\n", ind->obj_[obj_num - 1]);
            EP = CalChebycheff(ind, reference_point, obj_num);
            if (EP < bestEP)
            {
                bestEP = EP;
                bestIdx = i;
            }
            if (EP > worstEP)
            {
                worstEP = EP;
                worstIdx = i;
            }
        }
        fprintf(fpt3, "%f\n%f", bestEP, worstEP);


        Individual *ind_best = pop_table[bestIdx];
        Individual *ind_worst = pop_table[worstIdx];
        for (int j = 0; j < obj_num-1; ++j)
        {
            fprintf(fpt2, "%lf\t", reference_point[j]);
        }
        fprintf(fpt2, "%lf\n", reference_point[obj_num-1]);
        for (int j = 0; j < obj_num-1; ++j)
        {
            fprintf(fpt2, "%lf\t", ind_best->obj_[j]);
        }
        fprintf(fpt2, "%lf\n", ind_best->obj_[obj_num-1]);
        for (int j = 0; j < obj_num-1; ++j)
        {
            fprintf(fpt2, "%lf\t", ind_worst->obj_[j]);
        }
        fprintf(fpt2, "%lf\n", ind_worst->obj_[obj_num-1]);

        fclose(fpt1);
        fclose(fpt2);
        fclose(fpt3);
    }


    void RecordPop(int run_index, int generation, Global *para)
	{
        char output_dir[MAX_BUFFSIZE];
        char output_file_pop[MAX_BUFFSIZE];      // save whole population
        char para_settings[MAX_BUFFSIZE];
        // set the output directory
        std::string problem_name(para->problem_name_);
        std::string algorithm_name(para->algorithm_name_);
        // convert to uppercase
        for (auto &c : problem_name)
        {
            if (c >= '0' && c <= '9') continue;
            c = toupper(c);
        }
        for (auto &c : algorithm_name)
        {
            if (c >= '0' && c <= '9') continue;
            c = toupper(c);
        }

        if (para->algorithm_name_ == "moead_ltr" || para->algorithm_name_ == "r2_ibea_ltr")
        {
            sprintf(para_settings, "%d_%d_%.1f_%d_%.1f",
                    para->tau_,
                    para->inquiriesNum_,
                    para->step_size_,
                    MAX_QUEUE_SIZE,
                    para->kappa_
                    );
        }
        else
        {
            sprintf(para_settings, "%d_%d_%d_%.1f",
                    para->tau_,
                    para->inquiriesNum_,
                    MAX_QUEUE_SIZE,
                    para->kappa_
                    );
        }
        sprintf(output_dir, "../results");
        mkdir(output_dir,S_IRWXU);
        sprintf(output_dir, "../results/out");
        mkdir(output_dir,S_IRWXU);
        sprintf(output_dir, "../results/out/%s",
                algorithm_name.c_str());
        mkdir(output_dir,S_IRWXU);
        sprintf(output_dir, "../results/out/%s/%s",
                algorithm_name.c_str(),
                para_settings);
        mkdir(output_dir,S_IRWXU);
        sprintf(output_dir, "../results/out/%s/%s/%s_%dobjs_%ddvs",
                algorithm_name.c_str(),
                para_settings,
                para->problem_name_.c_str(),
                para->obj_num_,
                para->dec_num_);
        mkdir(output_dir,S_IRWXU);
        sprintf(output_dir, "../results/out/%s/%s/%s_%dobjs_%ddvs/%s",
                algorithm_name.c_str(),
                para_settings,
                para->problem_name_.c_str(),
                para->obj_num_,
                para->dec_num_,
                para->weight_stringType_.c_str());
        mkdir(output_dir,S_IRWXU);
        sprintf(output_dir, "../results/out/%s/%s/%s_%dobjs_%ddvs/%s/%d/",
                algorithm_name.c_str(),
                para_settings,
                para->problem_name_.c_str(),
                para->obj_num_,
                para->dec_num_,
                para->weight_stringType_.c_str(),
                run_index);
        mkdir(output_dir,S_IRWXU);

        sprintf(output_file_pop, "%spop_%d.txt", output_dir, generation);
        PrintObjective(output_file_pop, para->dec_num_, para->obj_num_, para->parent_population_.data(), para->population_num_, para->problem_name_);
	}

    void RecordPop(int run_index, int generation, Global *para, double *reference_point)
    {
        char output_dir_level0[MAX_BUFFSIZE];
        char output_dir_level1[MAX_BUFFSIZE];
        char output_dir_level2[MAX_BUFFSIZE];
        char output_dir_level3[MAX_BUFFSIZE];
        char output_dir_level4[MAX_BUFFSIZE];
        char output_file_pop[MAX_BUFFSIZE];      // save whole population
        char output_file_extremePop[MAX_BUFFSIZE];      // save extreme solutions
        char output_file_extremeEP[MAX_BUFFSIZE];      // save extreme EP

        // set the output directory

        std::string problem_name(para->problem_name_);
        std::string algorithm_name(para->algorithm_name_);
        for (auto &c : problem_name)
        {
            if (c >= '0' && c <= '9') continue;
            c = toupper(c);
        }
        for (auto &c : algorithm_name)
        {
            if (c >= '0' && c <= '9') continue;
            c = toupper(c);
        }

        
        sprintf(output_dir_level0, "../results/out");
        sprintf(output_dir_level1, "../results/out/%s",
                algorithm_name.c_str()
        );
        sprintf(output_dir_level2, "../results/out/%s/%s_%dobjs_%ddvs",
                algorithm_name.c_str(),
                para->problem_name_.c_str(),
                para->obj_num_,
                para->dec_num_
        );
        sprintf(output_dir_level3, "../results/out/%s/%s_%dobjs_%ddvs/%s/",
                algorithm_name.c_str(),
                para->problem_name_.c_str(),
                para->obj_num_,
                para->dec_num_,
                para->weight_stringType_.c_str()
        );
        sprintf(output_dir_level4, "../results/out/%s/%s_%dobjs_%ddvs/%s/%d/",
                algorithm_name.c_str(),
                para->problem_name_.c_str(),
                para->obj_num_,
                para->dec_num_,
                para->weight_stringType_.c_str(),
                run_index
        );
        #if defined(_WIN32)
                _mkdir(output_dir_level0);
                _mkdir(output_dir_level1);
                _mkdir(output_dir_level2);
        #elif defined(__linux) || defined(linux)
            mkdir(output_dir_level0,S_IRWXU);
            mkdir(output_dir_level1,S_IRWXU);
            mkdir(output_dir_level2,S_IRWXU);
            mkdir(output_dir_level3,S_IRWXU);
            mkdir(output_dir_level4,S_IRWXU);
        #endif
        sprintf(output_file_pop, "%spop_%d.txt", output_dir_level4, generation);
        sprintf(output_file_extremePop, "%sextremeSolution_%d.txt", output_dir_level4, generation);
        sprintf(output_file_extremeEP, "%sextremeEP_%d.txt", output_dir_level4, generation);
        PrintObjective(output_file_pop, output_file_extremePop, output_file_extremeEP, para->obj_num_, para->parent_population_.data(), para->population_num_, reference_point);
    }

	void ReadParametersFromFile(const char *filename, EMOCParameters *para)
	{
		int i = 0;
		char buff[MAX_BUFFSIZE] = { 0 };
		char line[MAX_BUFFSIZE] = { 0 };

		FILE *config = fopen(filename, "r");
		if (config == nullptr)
		{
			std::cerr << "Fail to open the config file, please check the file path." << std::endl;
			std::cerr << "Press enter to exit" << std::endl;
			std::cin.get();
			exit(-1);
		}

		while (!feof(config))
		{

			fgets(buff, MAX_BUFFSIZE, config);
			FormalizeStr(buff);
			for (i = 0; i < strlen(buff); i++)
			{
				if (buff[i] == ':')
				{
					buff[i] = 0;
					break;
				}
			}

			i++;
			if (!strcmp(buff, "algorithm_name"))
			{
				para->algorithm_name = buff + i;
			}
			else if (!strcmp(buff, "problem_name"))
			{
				para->problem_name = buff + i;
			}
			else if (!strcmp(buff, "variable_number"))
			{
				para->decision_num = atoi(buff + i);
			}
			else if (!strcmp(buff, "objective_number"))
			{
				para->objective_num = atoi(buff + i);
			}
			else if (!strcmp(buff, "population_number"))
			{
				para->population_num = atoi(buff + i);
			}
			else if (!strcmp(buff, "max_evaluation"))
			{
				para->max_evaluation = atoi(buff + i);
			}
			else if (!strcmp(buff, "output_interval"))
			{
				para->output_interval = atoi(buff + i);
			}
			else if (!strcmp(buff, "runs_number"))
			{
				para->runs_num = atoi(buff + i);
			}
			else if (!strcmp(buff, "open_multithread"))
			{
				para->is_open_multithread = atoi(buff + i);
			}
			else if (!strcmp(buff, "thread_number"))
			{
				para->thread_num = atoi(buff + i);
			}
            else if (!strcmp(buff, "askTime"))
            {
                para->askTime = atoi(buff + i);
            }
            else if (!strcmp(buff, "stepSize"))
            {
                para->step_size = atof(buff + i);
            }
            else if (!strcmp(buff, "elicitationInterval"))
            {
                para->tau = atoi(buff + i);
            }
            else if (!strcmp(buff, "inquiriesNum"))
            {
                para->inquiriesNum = atoi(buff + i);
            }
			else
			{
				std::cerr << "Input a wrong parameter, please check the parameter name" << std::endl;
				std::cerr << "Press enter to exit" << std::endl;
				std::cin.get();
				exit(-1);
			}

		}

		fclose(config);
	}

	void FormalizeStr(char *buff)
	{
		int i = 0, j = 0, len = 0;

		if (NULL == buff)
		{
			return;
		}

		len = strlen(buff);

		for (i = 0; i < len; i++)
		{
			switch (buff[i])
			{
			case ' ':
			case '\n':
				for (j = i; j < len; j++)
				{
					buff[j] = buff[j + 1];
				}
				break;
			case '\r':
				buff[i] = 0;
				break;
			default:
				break;
			}
		}
	}

	void ParseParamerters(int argc, char *argv[], EMOCParameters *para)
	{
		// defalut value
		para->algorithm_name = "nsga2";
		para->problem_name = "zdt1";
		para->population_num = 100;
		para->decision_num = 30;
		para->objective_num = 2;
		para->max_evaluation = 25000;
		para->output_interval = 100000;	// no output except the first and last gerneration
		para->runs_num = 1;
		para->is_open_multithread = 0;
		para->thread_num = 0;
		para->askTime = 10;
		para->tau = 0;
		para->inquiriesNum = 10;
		para->step_size = 0.5;
        para->kappa=-1;// no noise error

		// parse parameter from command line
		if (argc == 1)
			return;

		if (argc % 2 == 0)
		{
			std::cout << "The number of the parameter name is not matching the number of value." << std::endl;
			std::cout << "The parameters should be set like '-parameter_name value', e.g. '-algorithm nsga2'." << std::endl;
			std::cout << "Press enter to exit" << std::endl;
			std::cin.get();
			exit(-1);
		}
		else
		{
			char *first_parameter_name = argv[1] + 1;//skip '-'
            int i = 1;
			// read from file
			if (!strcmp(first_parameter_name, "file"))
			{
				ReadParametersFromFile(argv[2], para);
				i = 3;
			}
			//command line
            for (i; i+1 < argc;)
            {
                char *parameter_name = argv[i];//including '-'
                char *value = argv[i+1];
                SetParameter(parameter_name, value, para);
                i+=2;
            }
		}
	}
}