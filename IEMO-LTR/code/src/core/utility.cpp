#include "core/utility.h"

#include <cmath>

namespace emoc {

	DominateReleation CheckDominance(Individual *ind1, Individual *ind2, int obj_num)
	{
		int flag1 = 0, flag2 = 0;

		for (int i = 0; i < obj_num; ++i)
		{
			if (ind1->obj_[i] < ind2->obj_[i])
				flag1 = 1;
			else
			{
				if (ind1->obj_[i] > ind2->obj_[i])
					flag2 = 1;
			}
		}

		if (flag1 == 1 && flag2 == 0)
			return (DOMINATE);
		else
		{
			if (flag1 == 0 && flag2 == 1)
				return (DOMINATED);
			else
				return (NON_DOMINATED);
		}
	}

    DominateReleation check_r_dominance(Individual *ind1, Individual *ind2, int obj_num)
    {
        assert(obj_num == 2);

        double z[2] = {0};
        FILE *fpt = nullptr;
        fpt = fopen("/home/gylai/Desktop/ar_nsga2/code/w_approximate.txt", "r");
        if (!fpt)
        {
            std::cerr<<"<Error!!!> Could not open file"<<std::endl;
            exit(-1);
        }
        fscanf(fpt, "%lf %lf", &z[0], &z[1]);
        // std::cerr<<z[0]<<"\t"<<z[1]<<std::endl;
        fflush(fpt);
        fclose(fpt);

        double distance_x_z = 0;
        double distance_y_z = 0;

        int flag1 = 0, flag2 = 0;

        for (int i = 0; i < obj_num; ++i)
        {
            if (ind1->obj_[i] < ind2->obj_[i])
            {
                flag1 = 1;
            }
            else if(ind1->obj_[i] > ind2->obj_[i])
            {
                flag2 = 1;
            }

            distance_x_z += pow(ind1->obj_[i] - z[i], 2);
            distance_y_z += pow(ind2->obj_[i] - z[i], 2);
        }

        distance_x_z = sqrt(distance_x_z);
        distance_y_z = sqrt(distance_y_z);

        double r = .03;
        if ((flag1 == 1 && flag2 == 0)
            || (flag1 == 1 && flag2 == 1 && distance_x_z + r < distance_y_z))
        {
            return DOMINATE;
        }
        else if((flag1 == 0 && flag2 == 1)
                || (flag1 == 1 && flag2 == 1 && distance_x_z > distance_y_z + r))
        {
            return DOMINATED;
        }
        else
        {
            return NON_DOMINATED;
        }
    }


    int WeaklyDominates(double *point1, double *point2, int obj_num)
	{
		int i = 0, better = 1;
		while (i < obj_num && better)
		{
			better = point1[i] <= point2[i];
			i++;
		}
		return better;
	}

	double CalEuclidianDistance(double *a, double *b, int dimension)
	{
		double distance = 0.0;
		for (int i = 0; i < dimension; i++)
			distance += (a[i] - b[i]) * (a[i] - b[i]);
		return sqrt(distance);
	}

	double CalPerpendicularDistance(double *a, double *weight, int dimension)
	{
		double sin = CalculateSin(a, weight, dimension);
		double d2 = CalculateNorm(a, dimension);
		d2 = d2 * sin;

		return d2;
	}

	int Combination(int n, int k)
	{
		if (n < k)
			return -1;

		double ans = 1;
		for (int i = k + 1; i <= n; i++)
		{
			ans = ans * i;
			ans = ans / (double)(i - k);
		}

		return (int)ans;
	}

	double CalculateDotProduct(double *vector1, double *vector2, int dimension)
	{
		double dot_product = 0;
		for (int i = 0; i < dimension; i++)
			dot_product += vector1[i] * vector2[i];

		return dot_product;
	}

	double CalculateCos(double *a, double *b, int dimension)
	{
		return CalculateDotProduct(a, b, dimension) / (CalculateNorm(a, dimension) * CalculateNorm(b, dimension));
	}
    double CalculateCos(double *a, double *b, int dimension, int power)
    {
        return pow(CalculateDotProduct(a, b, dimension) / (CalculateNorm(a, dimension) * CalculateNorm(b, dimension)), power);
    }
	double L_alpha_norm(double *weight, double *obj, double *referencePoint, int dimension, double alpha)
	{
		double sum = 0.0;
		for (int i = 0; i < dimension; ++i)
		{
			//sum += pow(fabs(weight[i] * (obj[i] - referencePoint[i])), alpha);
			sum += pow(fabs((obj[i] - referencePoint[i]) / weight[i]), alpha);
		}
		return (pow(sum, 1.0 / alpha));
	}
	double L_alpha_norm_DM(double *weight, double *obj, int dimension, double alpha)
	{
		double sum = 0.0;
		for (int i = 0; i < dimension; ++i)
		{
			sum += pow(fabs(obj[i] / weight[i]), alpha);
		}
		return (pow(sum, 1.0 / alpha));
	}


	double CalculateSin(double *a, double *b, int dimension)
	{
		double cos = CalculateCos(a, b, dimension);
		return sqrt(1 - pow(cos, 2.0));
	}

	double CalculateNorm(double *vector, int dimension)
	{
		double norm = 0;
		for (int i = 0; i < dimension; i++)
		{
			norm += (vector[i] * vector[i]);
		}

		return sqrt(norm);
	}

	void UpdateIdealpoint(Individual *ind, double *ideal_point, int obj_num)
	{
		for (int i = 0; i < obj_num; i++)
		{
			if (ind->obj_[i] < ideal_point[i])
				ideal_point[i] = ind->obj_[i];
		}
	}

	void UpdateIdealpoint(Individual **pop, int pop_num, double *ideal_point, int obj_num)
	{
		for (int i = 0; i < obj_num; ++i)
			ideal_point[i] = INF;

		for (int i = 0; i < pop_num; i++)
		{
			for (int j = 0; j < obj_num; j++)
			{
				if (pop[i]->obj_[j] < ideal_point[j])
					ideal_point[j] = pop[i]->obj_[j];
			}
		}
	}

	void UpdateNadirpoint(Individual *ind, double *nadir_point, int obj_num)
	{
		for (int i = 0; i < obj_num; i++)
		{
			if (ind->obj_[i] > nadir_point[i])
				nadir_point[i] = ind->obj_[i];
		}
	}

	void UpdateNadirpoint(Individual **pop, int pop_num, double *nadir_point, int obj_num)
	{
		for (int i = 0; i < obj_num; ++i)
			nadir_point[i] = -INF;

		for (int i = 0; i < pop_num; i++)
		{
			for (int j = 0; j < obj_num; j++)
			{
				if (pop[i]->obj_[j] > nadir_point[j])
					nadir_point[j] = pop[i]->obj_[j];
			}
		}
	}

	double CalWeightedSum(Individual *ind, double *weight_vector, double *ideal_point, int obj_num)
	{
		double fitness = 0;
		for (int i = 0; i < obj_num; i++)
		{
			fitness += ind->obj_[i] * weight_vector[i];
		}

		ind->fitness_ = fitness;

		return fitness;
	}

	double CalInverseChebycheff(Individual *ind, double *weight_vector, double *ideal_point, int obj_num)
	{
		double fitness = 0, min = -1.0e+20;

		for (int i = 0; i < obj_num; ++i)
		{
			double diff = fabs(ind->obj_[i] - ideal_point[i]);
			if (weight_vector[i] < EPS)
				fitness = diff / 0.000001;
			else
				fitness = diff / weight_vector[i];

			if (fitness > min)
				min = fitness;
		}
		fitness = min;
		// ind->fitness_ = fitness;
		return fitness;
	}
    double CalInverseChebycheff(double *obj, double *weight_vector, double *ideal_point, int obj_num)
    {
        double fitness = 0, min = -1.0e+20;

        for (int i = 0; i < obj_num; ++i)
        {
            double diff = fabs(obj[i] - ideal_point[i]);
            if (weight_vector[i] < EPS)
                fitness = diff / 0.000001;
            else
                fitness = diff / weight_vector[i];

            if (fitness > min)
                min = fitness;
        }
        fitness = min;
        // ind->fitness_ = fitness;
        return fitness;
    }

    double CalInverseChebycheff(Individual *ind, double *weight_vector, int obj_num)
    {
        double fitness = 0, max = -1.0e+20;

        for (int i = 0; i < obj_num; ++i)
        {
            double diff = fabs(ind->obj_[i]);
            if (weight_vector[i] < EPS)
            {
                fitness = diff / 0.000001;
            }
            else
            {
                fitness = diff / weight_vector[i];
            }

            if (fitness > max)
            {
                max = fitness;
            }
        }
        fitness = max;
        return fitness;
    }

    double CalInverseChebycheff(double *obj, double *weight_vector, int obj_num)
    {
        double fitness = 0, max = -1.0e+20;

        for (int i = 0; i < obj_num; ++i)
        {
            double diff = fabs(obj[i]);
            if (weight_vector[i] < EPS)
            {
                fitness = diff / 0.000001;
            }
            else
            {
                fitness = diff / weight_vector[i];
            }

            if (fitness > max)
            {
                max = fitness;
            }
        }
        fitness = max;
        return fitness;
    }

    double CalChebycheff(Individual *ind, double *weight_vector, int obj_num)
    {
        double fitness = 0, max = -1.0e+20;

        for (int i = 0; i < obj_num; ++i)
        {
            double diff = fabs(ind->obj_[i]);
            if (weight_vector[i] < EPS)
            {
                fitness = diff * 0.000001;
            }
            else
            {
                fitness = diff * weight_vector[i];
            }

            if (fitness > max)
            {
                max = fitness;
            }
        }
        fitness = max;
        return fitness;
    }


	double CalPBI(Individual *ind, double *weight_vector, double *ideal_point, int obj_num, double theta)
	{
		theta == 0.0 ? 5.0 : theta;
		double d1 = 0.0, d2 = 0.0, nl = 0.0;

		for (int i = 0; i < obj_num; ++i)
		{
			d1 += (ind->obj_[i] - ideal_point[i]) * weight_vector[i];
			nl += weight_vector[i] * weight_vector[i];
		}
		nl = sqrt(nl);
		d1 = fabs(d1) / nl;

		for (int i = 0; i < obj_num; ++i)
			d2 += ((ind->obj_[i] - ideal_point[i]) - d1 * (weight_vector[i] / nl)) * ((ind->obj_[i] - ideal_point[i]) - d1 * (weight_vector[i] / nl));
		d2 = sqrt(d2);

		ind->fitness_ = d1 + theta * d2;
		return  (d1 + theta * d2);
	}

    Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd & origin, const float er = 0)
    {
        // 进行svd分解
        Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
                                                     Eigen::ComputeThinU |
                                                     Eigen::ComputeThinV);
        // 构建SVD分解结果
        Eigen::MatrixXd U = svd_holder.matrixU();
        Eigen::MatrixXd V = svd_holder.matrixV();
        Eigen::MatrixXd D = svd_holder.singularValues();

        // 构建S矩阵
        Eigen::MatrixXd S(V.cols(), U.cols());
        S.setZero();

        for (unsigned int i = 0; i < D.size(); ++i) {

            if (D(i, 0) > er) {
                S(i, i) = 1 / D(i, 0);
            }
            else {
                S(i, i) = 0;
            }
        }

        // pinv_matrix = V * S * U^T
        return V * S * U.transpose();
    }

/*calculate pinv | new zone | hold old zone on*/
    double** calculatePinv(double** n, int row, int col)
    {
        Eigen::MatrixXd originMatrix(row,col);
        Eigen::MatrixXd pinvMatrix;
        double** pinvN = nullptr;
        pinvN = new double*[col];
        for (int i = 0; i < row; ++i)
        {
            for (int j = 0; j < col; ++j)
            {
                originMatrix(i,j) = n[i][j];
            }
        }
        pinvMatrix = pinv_eigen_based(originMatrix);
        for (int i = 0; i < col; ++i)
        {
            pinvN[i] = new double[row];
            for (int j = 0; j < row; ++j)
            {
                pinvN[i][j] = pinvMatrix(i,j);
            }
        }
        return pinvN;
    }

    void TrainRBFNet(double **c, double *weight,double *out,  double sigma, int size, int dimension)
    {
        double **active_function = nullptr;
        double **active_function_pinv = nullptr;
        active_function = new double*[size];
        for(int i = 0; i < size; i++)
        {
            active_function[i] = new double[size];
        }
        double *center = nullptr;
        center = new double[dimension];
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                active_function[i][j] = 0;
            }
        }
        for (int i = 0;i < size;i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                center[j] = c[i][j];
            }
            for (int j = 0; j < size; j++)
            {
                for (int k = 0; k < dimension; k++)
                {
                    active_function[j][i] += pow(center[k] - c[j][k], 2);
                }
                active_function[j][i] = exp(sqrt(active_function[j][i])*(-sigma));
            }
        }
        active_function_pinv = calculatePinv(active_function, size, size);
        for (int i = 0; i < size; ++i)
        {
             weight[i] = 0.0;
            for (int j = 0; j < size; ++j)
            {
                weight[i] += out[j] * active_function_pinv[i][j];
            }
        }
        // memory free
        for (int i = 0; i < size; ++i)
        {
            delete[] active_function[i];
            delete[] active_function_pinv[i];
        }
        delete[] active_function;
        delete[] active_function_pinv;
        active_function = nullptr;
        active_function_pinv = nullptr;
        delete[] center;
        center = nullptr;
    }

    void UsingRBFNet(Individual *ind, double **c, double *weight, double sigma, int size, int dimension)
    {
        double *active_function, *center;
        active_function = center = nullptr;
        active_function = new double[size];
        center = new double[dimension];
        memset(active_function, 0, size * sizeof(double));
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < dimension; ++j)
            {
                center[j] = c[i][j];
            }
            for (int j = 0; j < dimension; ++j)
            {
                active_function[i] += pow(center[j] - ind->obj_[j], 2.0);
            }
            active_function[i] = exp(-1.0 * sqrt(active_function[i])*sigma);
        }
        ind->fitness_ = 0.0;
        for (int i = 0; i < size; ++i)
        {
            ind->fitness_ += weight[i] * active_function[i];
        }
        //free memory
        delete[] active_function;
        delete[] center;
        active_function = nullptr;
        center = nullptr;
    }

    double topK(double *model, double *gold, int len, int k)
    {
        /**
         * score via RankNet: the bigger, the better.
         */
        std::vector<double> scores_model(model,model+len);
        std::vector<size_t> idx_model(scores_model.size());
        iota(idx_model.begin(), idx_model.end(), 0);
        sort(idx_model.begin(), idx_model.end(),
             [&scores_model](size_t index_1, size_t index_2) { return scores_model[index_1] > scores_model[index_2]; });

        /**
         * score via CalInverseChebycheff: the smaller, the better.
         */
        std::vector<double> scores_gold(gold,gold+len);
        std::vector<size_t> idx_gold(scores_model.size());
        iota(idx_gold.begin(), idx_gold.end(), 0);
        sort(idx_gold.begin(), idx_gold.end(),
             [&scores_gold](size_t index_1, size_t index_2) { return scores_gold[index_1] < scores_gold[index_2]; });

        double res = 0.0;
        for (int i = 0; i < k; ++i)
        {
            for (int j = 0; j < k; ++j)
            {
                if (idx_model[i] == idx_gold[j])
                {
                    res = res + 1;
                    break;
                }
            }
        }
        res /= k;
        return res;
    }
    double* SetWeight(const std::string& weightstring)
    {
        int size = count(weightstring.begin(),weightstring.end(),',')+1;
        if(!size)
        {
            std::cout<<"\nthe dimensionality of weight is 0, hence exits."<<std::endl;
            exit(1);
        }
        double *m_w = nullptr;
        m_w = new double[size];
        int i = 0;
        //std::string::size_type -- unsigned long long
        std::string::size_type start = 0, end;
        while((end = weightstring.find(',',start)) != std::string::npos)
        {
            m_w[i++] = atof(weightstring.substr(start,end-start).c_str());
            start = end + 1;
        }
        m_w[i] = atof(weightstring.substr(start).c_str());
        return m_w;
    }

    void InitializeReferencePoint(double *reference_point, int obj_num)
    {
        for (int i = 0; i < obj_num; ++i)
        {
            reference_point[i] = 0;
        }
    }

    void UpdateReferencePoint(Individual **mixed_pop, int mixed_pop_num, int obj_num, double *reference_point)
    {
        double a,b, c;

        c = -INF;
        for (int i = 0; i < obj_num; ++i)
        {
            a = -INF;
            b = INF;
            for (int j = 0; j < mixed_pop_num; ++j)
            {
                if (mixed_pop[j]->obj_[i] > a)
                {
                    a = mixed_pop[j]->obj_[i];
                }
                if (mixed_pop[j]->obj_[i] < b)
                {
                    b = mixed_pop[j]->obj_[i];
                }
            }
            if ((a - b) > c)
            {
                c = a - b;
            }
        }
        #if DEBUG
        std::cout<<"reference_point: ";
        #endif
        for (int i = 0; i < obj_num; ++i)
        {
            a = INF;
            for (int j = 0; j < mixed_pop_num; ++j)
            {
                if (a > mixed_pop[j]->obj_[i])
                {
                    a = mixed_pop[j]->obj_[i];
                }
            }
            reference_point[i] = a - c;
            #if DEBUG
            std::cout<<reference_point[i]<<"\t";
            #endif
        }
        #if DEBUG
        std::cout<<"\n";
        #endif
    }

    void UpdateIdealPoint(Individual **pop, int pop_num, int obj_num, double *ideal_point)
    {
        for (int i = 0; i < pop_num; i++)
        {
            for (int j = 0; j < obj_num; j++)
            {
                if (pop[i]->obj_[j] < ideal_point[j])
                {
                    ideal_point[j] = pop[i]->obj_[j];
                }
            }
        }
        #if DEBUG
        std::cout<<"ideal_point: ";
        for (int i = 0; i < obj_num; ++i)
        {
            std::cout<<ideal_point[i]<<"\t";
        }
        std::cout<<"\n";
        #endif
    }

    float total_time = 0;
    void display_pop(FILE *gp,Individual **pop,int pop_num, int obj_num, int gen)
    {
        clock_t start, end;

        FILE *fptr= nullptr;
        fptr= fopen("plot.txt","w");
        if (!fptr)
        {
            std::cerr<<"<Error!!!> Could not open file"<<std::endl;
            exit(-1);
        }
        for (int i = 0; i < pop_num; i++)
        {
            for (int j = 0; j < obj_num-1; j++)
            {
                fprintf(fptr,"%f\t",pop[i]->obj_[j]);
            }
            fprintf(fptr,"%f\n",pop[i]->obj_[obj_num-1]);
        }
        fflush(fptr);
        fclose(fptr);
        char CMD[1024];
        if (obj_num==2)
        {
            sprintf(CMD,/*"set term post eps color enh solid\n"
                    "set term pdfcairo lw 2 font 'Times New Roman,8'\n"
                    "set output '../pop/pop_%d.pdf'\n"*/
                    "set grid\n"
                    "set autoscale\n"
                    "set title 'Generation #%d'\n"
                    "set xlabel 'f1'\n"
                    "set ylabel 'f2'\n"
                    "set xrange [0:1]\n"
                    "set yrange [-1:6]\n"
                    "unset key\n"
                    "plot 'PF_ZDT1.txt' w l lt -1 lw 2, 'plot.txt' w p pt 6 ps 1 lc 3, 'golden_point_zdt1_1_1.txt' w p pt 3 ps 2 lc 1\n"
                    ""
                    ,gen,gen);
        }
        else if (obj_num==3)
        {
            sprintf(CMD,"set term post eps color enh solid\n"
                    "set term pdfcairo lw 2 font 'Times New Roman,8'\n"
                    "set output '../pop/pop_%d.pdf'\n"
                    "set grid\n"
                    "set autoscale\n"
                    "set title 'Generation #%d'\n"
                    "set xlabel 'f1'\n"
                    "set ylabel 'f2'\n"
                    "set zlabel 'f3'\n"
                    "set xrange [0:1]\n"
                    "set yrange [0:1]\n"
                    "set zrange [0:1]\n"
                    "set view 45,45\n"
                    "unset key\n"
                    "splot 'plot.txt' w points pointtype 6 pointsize 1\n"
                    ,gen,gen);
        }
        else
        {
            std::cerr<<"Error!!! in display_pop(...)"<<std::endl;
            exit(-1);
        }


        start = clock();
        fprintf(gp,CMD);

        fflush(gp);

        end = clock();
        total_time += (double)(end - start)/CLOCKS_PER_SEC;
        std::cerr<<"run time:"<< (double)(end - start)/CLOCKS_PER_SEC << " total time:"<< total_time<<"\n";
        //show for a short time
        usleep(10000);
    }

    void display_mixed_pop(FILE *gp, int obj_num, int gen)
    {
        // FIXME: only support 2 objectives
        assert(obj_num == 2);
        char CMD[1024];
        sprintf(CMD,/*"set term post eps color enh solid\n"
                    "set term pdfcairo lw 2 font 'Times New Roman,8'\n"
                    "set output '../pop/pop_%d.pdf'\n"*/
                "set grid\n"
                "set autoscale\n"
                "set title 'Generation #%d'\n"
                "set xlabel 'f1'\n"
                "set ylabel 'f2'\n"
                "set xrange [0:1]\n"
                "set yrange [-1:6]\n"
                "unset key\n"
                "plot 'PF_ZDT1.txt' w l lt -1 lw 2,"
                "'mixed_pop.txt' w p pt 6 ps 1,"
                "'absolute_elite.txt' w p pt 7 ps 1,"
                "'survived_by_rank.txt' w p pt 9 ps 1,"
                "'golden_point_zdt1_1_1.txt' w p pt 3 ps 2\n"
                ,gen,gen);
        fprintf(gp,CMD);
        fflush(gp);
        //show for a short time
        usleep(10000);
    }

    void GNUPlot_score(FILE *gp,Individual **pop,int pop_num, int obj_num, int gen)
    {
//        FILE *fptr= nullptr;
//        fptr= fopen("plot.txt","w");
//        if (!fptr)
//        {
//            std::cerr<<"<Error!!!> Could not open file"<<std::endl;
//            exit(-1);
//        }
//        for (int i = 0; i < pop_num; i++)
//        {
//            for (int j = 0; j < obj_num-1; j++)
//            {
//                fprintf(fptr,"%f\t",pop[i]->obj_[j]);
//            }
//            fprintf(fptr,"%f\n",pop[i]->obj_[obj_num-1]);
//        }
//        fflush(fptr);
//        fclose(fptr);
        char CMD[1024];
        sprintf(CMD,"set term pdfcairo mono enhanced\n"
                    "set output '../GNUPlot_score/GNUPlot_score_%d.pdf'\n"
                    "set grid\n"
                    "set autoscale\n"
                    "set title 'Generation #%d'\n"
                    "set xlabel 'x'\n"
                    "set ylabel 'y'\n"
                    "set zlabel 'score'\n"
                    "set xrange [0:1]\n"
                    "set yrange [0:6]\n"
                    "set view 45,45\n"
                    "unset key\n"
                    "splot 'golden_function_1_1.txt' w p pt 7 ps 0.05\n"
                    "splot 'ranknet_score.txt' w p pt 1 ps 0.05\n"
                ,gen,gen);
        fprintf(gp,CMD);
        fflush(gp);
        //show for a short time
        //usleep(300000);
    }

    void initialize_bi_population(Individual **pop, int pop_num, int dec_num, int obj_num)
    {
        for (int i = 0; i < pop_num; ++i)
        {
            for (int j = 0; j < dec_num; ++j)
            {
                if (pop[i]->dec_[j] < 0.5)
                {
                    pop[i]->dec_[j] = 0;
                }else
                {
                    pop[i]->dec_[j] = 1;
                }
            }
        }
        greedy_repair(pop, pop_num, dec_num, obj_num);

    }

    bool isFeasible_knapsack(Individual *ind, int dec_num, int obj_num)
    {
        double total_weight = 0;
        for (int i = 0; i < obj_num; ++i)
        {
            total_weight = 0;
            for (int j = 0; j < dec_num; ++j)
            {
                total_weight += ind->dec_[j] * weight_knapsack[i][j];
            }
            if (total_weight > capacity_knapsack[i])
            {
                return false;
            }
        }
        return true;
    }

    void greedy_repair(Individual **pop, int pop_num, int dec_num, int obj_num)
    {
        for (int i = 0; i < pop_num; ++i)
        {
            // greedy repair method
            greedy_repair_ind(pop[i], dec_num, obj_num);
        }
    }

    void greedy_repair_ind(Individual *ind, int dec_num, int obj_num)
    {
        // check the constraints for knapsack and repair
        while (!isFeasible_knapsack(ind, dec_num, obj_num))
        {
            double min_temp = INF;
            int idx_min = 0;
            for (int i = 0; i < dec_num; ++i)
            {
                if (fabs(ind->dec_[i] - 1) < 1e-6)
                {// check for one item which has been selected
                    double max_temp = -1;
                    double temp = 0;
                    for (int j = 0; j < obj_num; ++j)
                    {
                        temp = static_cast<double>(profit_knapsack[j][i]) / static_cast<double>(weight_knapsack[j][i]);
                        if(temp > max_temp)
                        {
                            max_temp = temp;
                        }
                    }
                    if (max_temp < min_temp)
                    {
                        min_temp = max_temp;
                        idx_min = i;
                    }
                }
            }
            ind->dec_[idx_min] = 0;
        }
    }

    int get_random_integer(int min, int max)
    {
        return (rand() % (max-min+1))+ min;
    }

}