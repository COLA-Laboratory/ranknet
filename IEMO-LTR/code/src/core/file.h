#include "core/global.h"
#include "core/individual.h"

namespace emoc {

	struct EMOCParameters
	{
		std::string algorithm_name;
		std::string problem_name;
		int population_num;
		int decision_num;
		int objective_num;
		int max_evaluation;
		int output_interval;
		int start_index;
		int runs_num;
		int is_open_multithread;
		int thread_num;
		int askTime;    //total time of consultation
		int tau; // elicitation interval
		int inquiriesNum;   //number of pairwise comparisons in one consultation session
        std::string weight_stringType;
        double kappa;
        double step_size;
		double *igd_value;
	};

	void PrintObjective(const char *filename, int obj_num, Individual **pop_table, int pop_num, const std::string &problem_name);
	void PrintObjective(const char *filename_pop, const char *filename_extremeSolution, const char *filename_extremeEP, int obj_num, Individual **pop_table, int pop_num, double *goldenPoint);
	void RecordPop(int run_index, int generation, Global *para);
    void RecordPop(int run_index, int generation, Global *para, double *reference_point);
	void ReadParametersFromFile(const char *filename, EMOCParameters *para);
	void ParseParamerters(int argc, char *argv[], EMOCParameters *para);
	void FormalizeStr(char *buff);
}