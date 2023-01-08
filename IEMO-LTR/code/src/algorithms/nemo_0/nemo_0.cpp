//
// Created by guiyulai on 2021/5/17.
//
#include "algorithms/nemo_0/nemo_0.h"

namespace emoc {

    NEMO_0::NEMO_0(Problem *problem, int thread_num):
    Algorithm(problem, thread_num),
    weight_(nullptr)
    {

    }

    NEMO_0::~NEMO_0()
    {
        if (weight_)
        {
            delete[] weight_;
            weight_ = nullptr;
        }
        while(rank_of_winners_.size())
        {
            delete[] rank_of_winners_.front();
            rank_of_winners_.pop();
            delete[] rank_of_losers_.front();
            rank_of_losers_.pop();
        }
    }

    void NEMO_0::Run()
    {

        STDCOUT("NEMO_0::Run()");

        int tau;
        tau = g_GlobalSettings->tau_;// TODO
        bool flag_EnvironmentalSelection;
        Initialization();

        if (g_GlobalSettings->output_interval_ != 0)
        {
            TrackPopulation(g_GlobalSettings->iteration_num_);
        }
        while (!g_GlobalSettings->IsTermination())
        {
            STDCOUT("begin each iteration");
            // begin each iteration
            g_GlobalSettings->iteration_num_++;
            // generate offspring population
            Crossover(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->offspring_population_.data());
            MutationPop(g_GlobalSettings->offspring_population_.data(), 2 * g_GlobalSettings->population_num_ / 2, g_GlobalSettings);
            EvaluatePop(g_GlobalSettings->offspring_population_.data(), 2 * g_GlobalSettings->population_num_ / 2);
            MergePopulation(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_, g_GlobalSettings->offspring_population_.data(),
                            2 * g_GlobalSettings->population_num_ / 2, g_GlobalSettings->mixed_population_.data());

            STDCOUT("ready to EnvironmentalSelection");
            flag_EnvironmentalSelection = EnvironmentalSelection(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->mixed_population_.data(), (g_GlobalSettings->iteration_num_ % tau == 0));
            if (!flag_EnvironmentalSelection)
            {
                return;
            }
            // record the population every interval generations and the first and last genration
            if (((g_GlobalSettings->output_interval_ != 0)&&
            (g_GlobalSettings->iteration_num_ % g_GlobalSettings->output_interval_ == 0))
            || g_GlobalSettings->IsTermination())
            {
                TrackPopulation(g_GlobalSettings->iteration_num_);
            }
        }
    }

    void NEMO_0::Initialization()
    {
        // initialize parent population
        g_GlobalSettings->InitializePopulation(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);
        EvaluatePop(g_GlobalSettings->parent_population_.data(), g_GlobalSettings->population_num_);
        // set weight
        weight_ = SetWeight(g_GlobalSettings->weight_stringType_);
    }

    void NEMO_0::Crossover(Individual **parent_pop, Individual **offspring_pop)
    {
        // generate random permutation index for tournment selection
        int *index1 = new int[g_GlobalSettings->population_num_];
        int *index2 = new int[g_GlobalSettings->population_num_];
        random_permutation(index1, g_GlobalSettings->population_num_);
        random_permutation(index2, g_GlobalSettings->population_num_);

        for (int i = 0; i < g_GlobalSettings->population_num_ / 2; ++i)
        {
            Individual *parent1 = TournamentByRank(parent_pop[index1[2 * i]], parent_pop[index1[2 * i + 1]]);
            Individual *parent2 = TournamentByRank(parent_pop[index2[2 * i]], parent_pop[index2[2 * i + 1]]);
            SBX(parent1, parent2, offspring_pop[2 * i], offspring_pop[2 * i + 1], g_GlobalSettings);
        }
        if (index1)
        {
            delete[] index1;
            index1 = nullptr;
        }
        if (index2)
        {
            delete[] index2;
            index2 = nullptr;
        }
    }

    void NEMO_0::SetDistanceInfo(std::vector<DistanceInfo> &distanceinfo_vec, int target_index, double distance)
    {
        // search the target_index and set it's distance
        for (int i = 0; i < distanceinfo_vec.size(); ++i)
        {
            if (distanceinfo_vec[i].index == target_index)
            {
                distanceinfo_vec[i].distance += distance;
                break;
            }
        }
    }
    void NEMO_0::SetPreferenceInfo(std::vector<SortList> &preferenceInfo_vec, int target_index, double value)
    {
        // search the target_index and set it's distance
        for (int i = 0; i < preferenceInfo_vec.size(); ++i)
        {
            if (preferenceInfo_vec[i].index == target_index)
            {
                preferenceInfo_vec[i].value += value;
                break;
            }
        }
    }


    int NEMO_0::CrowdingDistance(Individual **mixed_pop,  int pop_num, int *pop_sort, int rank_index)
    {
        int num_in_rank = 0;
        std::vector<int> sort_arr(pop_num, 0);
        std::vector<DistanceInfo> distanceinfo_vec(pop_num, { -1,0.0 });

        // find all the indviduals with rank rank_index
        for (int i = 0; i < pop_num; i++)
        {
            mixed_pop[i]->fitness_ = 0;
            if (mixed_pop[i]->rank_ == rank_index)
            {
                distanceinfo_vec[num_in_rank].index = i;
                sort_arr[num_in_rank] = i;
                num_in_rank++;
            }
        }
        for (int i = 0; i < g_GlobalSettings->obj_num_; i++)
        {
            // sort the population with i-th obj
            std::sort(sort_arr.begin(), sort_arr.begin()+num_in_rank, [=](const int left, const int right){
                return mixed_pop[left]->obj_[i] < mixed_pop[right]->obj_[i];
            });

            // set the first and last individual with INF fitness (crowding distance)
            mixed_pop[sort_arr[0]]->fitness_ = INF;
            SetDistanceInfo(distanceinfo_vec, sort_arr[0], INF);
            mixed_pop[sort_arr[num_in_rank - 1]]->fitness_ = INF;
            SetDistanceInfo(distanceinfo_vec, sort_arr[num_in_rank - 1], INF);

            // calculate each solution's crowding distance
            for (int j = 1; j < num_in_rank - 1; j++)
            {
                if (INF != mixed_pop[sort_arr[j]]->fitness_)
                {
                    if (mixed_pop[sort_arr[num_in_rank - 1]]->obj_[i] == mixed_pop[sort_arr[0]]->obj_[i])
                    {
                        mixed_pop[sort_arr[j]]->fitness_ += 0;
                    }
                    else
                    {
                        double distance = (mixed_pop[sort_arr[j + 1]]->obj_[i] - mixed_pop[sort_arr[j - 1]]->obj_[i]) /
                                          (mixed_pop[sort_arr[num_in_rank - 1]]->obj_[i] - mixed_pop[sort_arr[0]]->obj_[i]);
                        mixed_pop[sort_arr[j]]->fitness_ += distance;
                        SetDistanceInfo(distanceinfo_vec, sort_arr[j], distance);
                    }
                }
            }
        }


        std::sort(distanceinfo_vec.begin(), distanceinfo_vec.begin()+num_in_rank, [](const DistanceInfo &left, const DistanceInfo &right) {
            return left.distance < right.distance;
        });

        // copy sort result
        for (int i = 0; i < num_in_rank; i++)
        {
            pop_sort[i] = distanceinfo_vec[i].index;
        }

        return num_in_rank;
    }
    int NEMO_0::ArbitrarilySelect(Individual **mixed_pop, int pop_num, int *pop_sort, int rank_index)
    {
        int num_in_rank = 0;
        std::vector<int> sort_arr(pop_num, 0);
        std::vector<DistanceInfo> distanceinfo_vec(pop_num, { -1,0.0 });

        // find all the indviduals with rank rank_index
        for (int i = 0; i < pop_num; i++)
        {
            mixed_pop[i]->fitness_ = 0;
            if (mixed_pop[i]->rank_ == rank_index)
            {
                distanceinfo_vec[num_in_rank].index = i;
                sort_arr[num_in_rank] = i;
                num_in_rank++;
            }
        }
        std::sort(distanceinfo_vec.begin(), distanceinfo_vec.begin()+num_in_rank, [](const DistanceInfo &left, const DistanceInfo &right) {
            return left.distance < right.distance;
        });

        // copy sort result
        for (int i = 0; i < num_in_rank; i++)
        {
            pop_sort[i] = distanceinfo_vec[i].index;
        }
        return num_in_rank;
    }

    int NEMO_0::RankViaPreference(Individual **mixed_pop, int pop_num, int *pop_sort, int rank_index)
    {
        STDCOUT("NEMO_0::RankViaPreference");
        int number_of_inquiries, index_1, index_2, index_of_ind_1, index_of_ind_2, idx2ind;
        int *rank_1, *rank_2;
        bool isSuccess = false;
        rank_1 = rank_2 = nullptr;
        int size;
        int num_in_rank = 0;
        std::vector<int> PF(pop_num, 0);// pop_num: 2 * pop_parent_num
        std::vector<int> lp_rank(pop_num, 0);
        std::vector<SortList> preferenceInfo_vec(pop_num, { -1,0.0 });

        // find all the indviduals with rank rank_index
        for (int i = 0; i < pop_num; i++)
        {
            if (mixed_pop[i]->rank_ == rank_index)
            {
                preferenceInfo_vec[num_in_rank].index = i;
                PF[num_in_rank] = i;
                num_in_rank++;
            }
        }
        size = num_in_rank;
        number_of_inquiries = g_GlobalSettings->inquiriesNum_;
        while (number_of_inquiries--)
        {
            if (rank_of_winners_.size()==MAX_QUEUE_SIZE)
            {
                rank_1 = rank_of_winners_.front();
                if (rank_1)
                {
                    delete[] rank_1;
                    rank_1 = nullptr;
                }
                rank_of_winners_.pop();
                rank_2 = rank_of_losers_.front();
                if (rank_2)
                {
                    delete[] rank_2;
                    rank_2 = nullptr;
                }
                rank_of_losers_.pop();
            }

            index_1 = index_2 = 0;
            while((index_1==index_2))
            {
                index_1 = rnd(0, num_in_rank - 1);
                index_2 = rnd(0, num_in_rank - 1);
            }
            index_of_ind_1 = PF[index_1];
            index_of_ind_2 = PF[index_2];
            // calculate the rank of objectives in every dimension for this two individuals
            rank_1 = new int[g_GlobalSettings->obj_num_];
            rank_2 = new int[g_GlobalSettings->obj_num_];
            for (int i = 0; i < g_GlobalSettings->obj_num_; ++i)
            {
                rank_1[i] = rank_2[i] = 0;
                for (int j = 0; j < num_in_rank; ++j)
                {
                    idx2ind = PF[j];
                    if (g_GlobalSettings->mixed_population_[index_of_ind_1]->obj_[i]
                        >= g_GlobalSettings->mixed_population_[idx2ind]->obj_[i])
                    {
                        rank_1[i]++;
                    }
                    if (g_GlobalSettings->mixed_population_[index_of_ind_2]->obj_[i]
                        >= g_GlobalSettings->mixed_population_[idx2ind]->obj_[i])
                    {
                        rank_2[i]++;
                    }
                }
            }

            if (CalInverseChebycheff(g_GlobalSettings->mixed_population_[index_of_ind_1], weight_, g_GlobalSettings->obj_num_)
                <CalInverseChebycheff(g_GlobalSettings->mixed_population_[index_of_ind_2], weight_, g_GlobalSettings->obj_num_))
            {
                rank_of_winners_.push(rank_1);
                rank_of_losers_.push(rank_2);
            }
            else
            {
                rank_of_winners_.push(rank_2);
                rank_of_losers_.push(rank_1);
            }
        }
        // solving LP
        double avf[num_in_rank * g_GlobalSettings->obj_num_ + 2];
        isSuccess = LP_Solver(avf, size, g_GlobalSettings->obj_num_);

        if (isSuccess)
        {
            #if DEBUG
            std::cout<<"[INFO] Solve lp successfully."<<std::endl;
            #endif
        }
        else
        {
            #if DEBUG
            std::cout<<"[ERROR] "<<g_GlobalSettings->problem_name_<<"\t"<<g_GlobalSettings->weight_stringType_<<std::endl;
            std::cout<<"[ERROR] Solve lp unsuccessfully in run No."<<g_GlobalSettings->run_id_<<std::endl;
            #endif
            std::cerr<<"[ERROR] "<<g_GlobalSettings->problem_name_<<"\t"<<g_GlobalSettings->weight_stringType_<<std::endl;
            std::cerr<<"[ERROR] Solve lp unsuccessfully in run No."<<g_GlobalSettings->run_id_<<std::endl;
            return -1;
        }
        #if DEBUG
        std::cout<<"[INFO] The objective values of the last PF:\n";
        for (int i = 0; i < num_in_rank; ++i)
        {
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                std::cout<<g_GlobalSettings->mixed_population_[PF[i]]->obj_[j]<<"\t";
            }
            std::cout<<std::endl;
        }
        #endif
        // get rank of each objective in every dimension for the last PF
        int ranking[num_in_rank][g_GlobalSettings->obj_num_];
        for (int i = 0; i < num_in_rank; ++i)
        {
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                ranking[i][j] = 0;
                for (int k = 0; k < num_in_rank; ++k)
                {
                    if (g_GlobalSettings->mixed_population_[PF[i]]->obj_[j]
                        >= g_GlobalSettings->mixed_population_[PF[k]]->obj_[j])
                    {
                        ranking[i][j]++;
                    }
                }
            }
        }
        /*for (int i = 0; i < g_GlobalSettings->obj_num_; ++i)
        {
            for (int j = 0; j < num_in_rank; ++j)
            {
                lp_rank[j] = j + 1;
            }
            std::sort(lp_rank.begin(), lp_rank.end(),
                  [=](const int left, const int right)
                    {
                        return mixed_pop[PF[left]]->obj_[i] < mixed_pop[PF[right]]->obj_[i];
                    });
            for (int j = 0; j < num_in_rank; ++j)
            {
                ranking[j][i] = lp_rank[j];// TODO: DEBUG
            }
        }*/
        #if DEBUG
        std::cout<<"[INFO] Rank of the last PF in Gen NO."<<g_GlobalSettings->iteration_num_<<":"<<std::endl;
        for (int i = 0; i < num_in_rank; ++i)
        {
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                std::cout<<ranking[i][j]<<"\t";
            }
            std::cout<<std::endl;
        }
        #endif
        int idx;
        int idx_ind;
        for (int i = 0; i < num_in_rank; ++i)
        {
            idx_ind = PF[i];
            mixed_pop[idx_ind]->fitness_ = 0.0;
            for (int j = 0; j < g_GlobalSettings->obj_num_; ++j)
            {
                idx = j * size + ranking[i][j] - 1;
                mixed_pop[idx_ind]->fitness_ += avf[idx];
                SetPreferenceInfo(preferenceInfo_vec, idx_ind, avf[idx]);
            }
        }
        #if DEBUG
        std::cout<<"[INFO] Score in the last PF:\n";
        for (int i = 0; i < num_in_rank; ++i)
        {
            std::cout<<preferenceInfo_vec[i].index<<": "<<preferenceInfo_vec[i].value<<std::endl;
        }
        #endif
        #if UsingGoldenFunction
        std::cout<<"[INFO] Using golden function.\n";
        for (int i = 0; i < num_in_rank; ++i)
        {
            preferenceInfo_vec[i].value = 0.0;
        }
        for (int i = 0; i < num_in_rank; ++i)
        {
            idx_ind = PF[i];
            SetPreferenceInfo(preferenceInfo_vec, idx_ind, 1-CalInverseChebycheff(mixed_pop[idx_ind], weight_, g_GlobalSettings->obj_num_));
        }
        #endif
        std::sort(preferenceInfo_vec.begin(), preferenceInfo_vec.begin()+num_in_rank, [](const SortList &left, const SortList &right) {
            return left.value < right.value;
        });
        #if DEBUG
        std::cout<<"[INFO] Rank of the last PF:\n";
        for (int i = 0; i < num_in_rank; ++i)
        {
            std::cout<<preferenceInfo_vec[i].index<<": "<<preferenceInfo_vec[i].value<<std::endl;
        }
        #endif

        // copy sort result
        for (int i = 0; i < num_in_rank; i++)
        {
            pop_sort[i] = preferenceInfo_vec[i].index;
        }
        return num_in_rank;
    }

    bool NEMO_0::EnvironmentalSelection(Individual **parent_pop, Individual **mixed_pop, bool usingPreferenceModel)
    {
        STDCOUT("NEMO_0::EnvironmentalSelection");
        int current_popnum = 0, rank_index = 0;
        int mixed_popnum = g_GlobalSettings->population_num_ + 2 * g_GlobalSettings->population_num_ / 2;//TODO: ???
        // do nondominated sorting in the mixed population as the first selection pressure
        NonDominatedSort(mixed_pop, mixed_popnum, g_GlobalSettings->obj_num_);
        // select individuals by rank
        while (1)
        {
            int temp_number = 0;// record the number of Pareto front of rank_index
            for (int i = 0; i < mixed_popnum; i++)
            {
                if (mixed_pop[i]->rank_ == rank_index)
                {
                    temp_number++;
                }
            }
            if (current_popnum + temp_number <= g_GlobalSettings->population_num_)
            {// add all individuals in Pareto front of rank_index to the population of the new offspring
                for (int i = 0; i < mixed_popnum; i++)
                {
                    if (mixed_pop[i]->rank_ == rank_index)
                    {
                        CopyIndividual(mixed_pop[i], parent_pop[current_popnum]);
                        current_popnum++;
                    }
                }
                rank_index++;
            }
            else// only rely on the 1st selection pressure is not enough, so elicit the crowding distance as the second selection pressure
            {
                break;
            }
        }

        // select individuals by crowding distance
        int sort_num = 0;// record
        int *pop_sort = new int[mixed_popnum];

        if (current_popnum < g_GlobalSettings->population_num_)
        {
            if (!usingPreferenceModel)
            {
                sort_num = ArbitrarilySelect(mixed_pop, mixed_popnum, pop_sort, rank_index);
            }
            else
            {
                sort_num = RankViaPreference(mixed_pop, mixed_popnum, pop_sort, rank_index);
                if (sort_num == -1)
                {
                    return false;
                }
            }
            while (1)
            {
                if (current_popnum < g_GlobalSettings->population_num_)
                {
                    CopyIndividual(mixed_pop[pop_sort[--sort_num]], parent_pop[current_popnum]);//(A, B): A->B
                    current_popnum++;
                }
                else
                {
                    break;
                }
            }
        }

        // clear crowding distance value
        for (int i = 0; i < g_GlobalSettings->population_num_; i++)
        {
            parent_pop[i]->fitness_ = 0.0;
        }
        if (pop_sort)
        {
            delete[] pop_sort;
            pop_sort = nullptr;
        }
        return true;
    }

    bool NEMO_0::LP_Solver(double * row, int size, int obj_num)
    {
        STDCOUT("NEMO_0::LP_Solver");

        #if DEBUG
        std::cout<<"[INFO] size: "<<size<<"\tobj_num: "<<obj_num<<std::endl;
        #endif
        int pairs_num = rank_of_winners_.size();
        int variable_num = size * obj_num + 1;
        #if DEBUG
        std::cout<<"[INFO] pairs_num: "<<pairs_num<<"\tvariable_num: "<<variable_num<<std::endl;
        #endif
        int constraints_num = 0;
        int ret = 0;
        int t;
        int winner_index, loser_index;
        int winner[obj_num], loser[obj_num];
        //lp_solve model
        lprec *lp;
        /* must be 1 more than number of columns ! */
        //create new LP model row by row with 0 rows and number_variable columns
        lp = make_lp(0, variable_num);
        if(lp == NULL)
        {
            std::cerr<<"[ERROR] Unable to create new LP model."<<std::endl;
            return false;
        }
        else
        {
            #if DEBUG
            std::cout<<"[INFO] Create new LP model successfully."<<std::endl;
            #endif
        }
        //add constraints row by row
        set_add_rowmode(lp, TRUE);
        //set the object direction to maximize
        set_maxim(lp);
        #if DEBUG
        std::cout<<"[INFO] Add constraints to make sure that for each objective u_i, u_i(rank smaller) >= u_i(rank_bigger)."<<std::endl;
        #endif
        //add constraint
        for (int i = 0; i < obj_num; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                t = size - j - 1;
                for (int k = 0; k < t; ++k)
                {
                    //row --> double matrix, size: 1 * (1 + number_variable)
                    for (int m = 1; m < variable_num + 1; m++)
                    {
                        row[m] = 0;
                    }
                    row[size * i + j +1] = -1;
                    row[size * i + j + k + 2] = 1;
                    /* reference: http://lpsolve.sourceforge.net/5.5/
                     *
                     * Add a constraint to the lp.
                     * unsigned char add_constraint(lprec *lp, REAL *row, int constr_type, REAL rh);
                     * lp: Pointer to previously created lp model.
                     * row: An array with 1+get_Ncolumns (count for add_constraintex, if colno is different from NULL) elements that contains the values of the row.
                     * constr_type: LE (1) --> <=
                     *              EQ (3) --> =
                     *              GE (2) --> >=
                     * rh: The value of the right hand side (RHS).
                     * return: TRUE (1) if the operation was successful. A return value of FALSE (0) indicates an error.
                     *
                     * Remarks: Note that for add_constraint (and add_constraintex when colno is NULL)
                     * element 0 of the array is not considered (i.e. ignored).
                     * Column 1 is element 1, column 2 is element 2, ...
                     * str_add_constraint should only be used in small or demo code
                     * since it is not performant and uses more memory.
                     *
                     * */
                    add_constraint(lp, row, LE, 0);//LE --> 1
                }
                constraints_num += t;
            }
        }
        /**
         * x_win - x_loss >= epsilon
         */
        for (int i = 0; i < pairs_num; ++i)
        {
            #if DEBUG
            std::cout<<"[INFO] Pair NO."<<i<<std::endl;
            #endif
            for (int j = 0; j < obj_num; ++j)
            {
                winner[j] = (rank_of_winners_.front())[j];
                #if DEBUG
                std::cout<<winner[j]<<"\t";
                #endif
            }
            rank_of_winners_.push(rank_of_winners_.front());
            rank_of_winners_.pop();
            #if DEBUG
            std::cout<<std::endl;
            #endif
            for (int j = 0; j < obj_num; ++j)
            {
                loser[j] = (rank_of_losers_.front())[j];
                #if DEBUG
                std::cout<<loser[j]<<"\t";
                #endif
            }
            std::cout<<std::endl;
            rank_of_losers_.push(rank_of_losers_.front());
            rank_of_losers_.pop();
            for (int m = 1; m < variable_num + 1; m++)
            {
                row[m] = 0;
            }
            row[variable_num]=1;
            for (int j = 0;j < obj_num;j++)
            {
                winner_index = winner[j] + size * j;
                row[winner_index] = -1;
                loser_index = loser[j] + size * j;
                row[loser_index] = 1;
            }
            add_constraint(lp, row, LE, 0);
            constraints_num++;
        }
        /* construct the objective function for lp: max epsilon */
        for (int m = 1; m < variable_num + 1; m++)
        {
            row[m] = 0;//ignore row[0]
        }
        row[variable_num] = 1;//number_variable --> front_size * number_objective + 1
        /* reference: http://lpsolve.sourceforge.net/5.5/
         *
         * unsigned char set_obj_fn(lprec *lp, REAL *row);
         * Set the objective function (row 0) of the matrix.
         * return TRUE (1) if the operation was successful.
         * A return value of FALSE (0) indicates an error.
         * row: An array with 1+get_Ncolumns (count for set_obj_fnex)
         * elements that contains the values of the objective function.
         * The set_obj_fn functions set all values of the objective
         * function at once.
         * Note that for set_obj_fn element 0 of the array is not
         * considered (i.e. ignored). Column 1 is element 1, column 2 is element 2, ...
         *
         * Note that it is advised to set the objective function
         * before adding rows via add_constraint, add_constraintex,
         * str_add_constraint. This especially for larger models.
         * This will be much more performant than adding the objective
         * function afterwards.
         * */
        set_obj_fn(lp, row);//ignore row[0]
        /**
         * u_1(worst) = u_2(worst) = ... = u_m(worst) = 0
         */
        for (int m = 1; m < variable_num + 1; m++)
        {
            row[m] = 0;
        }
        for (int i = 0;i < obj_num;i++)
        {
            row [size*(i+1)] = 1;
        }
        add_constraint(lp, row, EQ, 0);
        constraints_num++;
        /**
         * u_1(best) + u_2(best) + ... + u_m(best) = 1
         */
        for (int m = 1; m < variable_num + 1; m++)
        {
            row[m] = 0;
        }
        for (int i = 0;i < obj_num;i++)
        {
            row [size * i + 1 ] = 1;
        }
        add_constraint(lp, row, EQ, 1);
        constraints_num++;
        for (int i = 1 ; i < variable_num + 1;i++)
        {
            set_bounds(lp, i, 0.0,1.0);
        }
        #if DEBUG
        std::cout<<"[INFO] Number of all constraints in LP: "<<constraints_num<<std::endl;
        #endif
        /* rowmode should be turned off again when done building the model */
        set_add_rowmode(lp, FALSE);//close the switch
        /* I only want to see important messages on screen while solving */
        set_verbose(lp, IMPORTANT);
        /* Now let lpsolve calculate a solution */
        ret = solve(lp);//solve the LP

        if (ret == OPTIMAL)//OPTIMAL --> 0
        {
            #if DEBUG
            std::cout<<"[INFO] Obtained the optimal solution in lp_solver"<<std::endl;
            #endif
            ret = 0;
        }
        else
        {
            std::cerr<<"[ERROR] Lp_solver did not obtain the optimal solution"<<std::endl;
            ret = 5;
            return false;
        }
        if (ret == 0)/* a solution is calculated, now let's get some results */
        {
            /* objective value */
            #if DEBUG
            std::cout<<"[INFO] Objective value: "<<get_objective(lp)<<std::endl;
            #endif
            /* variable values */
            /* reference: http://lpsolve.sourceforge.net/5.5/
             * unsigned char get_variables(lprec *lp, REAL *var);
             * Returns the values of the variables.
             * lp: Pointer to previously created lp model.
             * var: An array that will contain the values of the variables.
             * returns TRUE (1) if the operation was successful. A return value of FALSE (0) indicates an error.
             * Note: Element 0 will contain the value of the first variable, element 1 of the second variable, ...
             * */
            get_variables(lp, row);//get the result, saving in the matrix $row
            #if DEBUG
            for(int j = 0; j < variable_num; j++)
            {
                //row[number_variable-1] is the objective value of solved lp
                std::cout<<"[INFO] "<<get_col_name(lp, j + 1)<<": "<<row[j]<<std::endl;
            }
            #endif
        }
        if (lp !=nullptr)
        {
            /* clean up such that all used memory by lpsolve is freed */
            delete_lp(lp);
        }

        STDCOUT("NEMO_0::LP_Solver returns.");

        return true;
    }


}

