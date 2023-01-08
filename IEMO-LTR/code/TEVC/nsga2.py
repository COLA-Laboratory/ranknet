#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Public.InitialPop import initial
from Non_Dominated_Sort import nonDominatedSort
from CrowdingDistance import crowdingDistance
from NSGAIISelection import tournamentSelection, elitism
from Operator.CrossAndMutation import crossMutation
from ProblemSelect import problemSelect
from latent_preference import golden_function

# RankNet related
from RankNet import CreateModel
from RankNet import TrainModel
from RankNet import UseModel

# In[2]:


import time
import numpy as np
import random
import csv
import math
import matplotlib.pyplot as plt

# # parameters

# In[3]:

def main():
    problemName = 'swimmer'
    N = 100
    maxgen = 300
    D = 18
    M = 2

    # # NSGA-II

    # In[4]:


    problem = problemSelect(problemName, D, M)
    start = time.time()

    D = problem.D
    M = problem.M
    lower = problem.lower
    upper = problem.upper
    pc = 0.9
    pm = 1 / D
    pop = initial(N, D, M, lower, upper, problem, 'real')

    F1, pop_non, _ = nonDominatedSort(N, pop)
    pop = crowdingDistance(F1, pop_non, M)
    gen = 1

    # preference related
    tau = 10
    mu = 10
    winners = []
    losers = []


    # # results

    # In[6]:

    systemTime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    f = open("logs/log%f.csv" % time.time(), 'a', encoding='utf-8', newline='')
    csv_writer = csv.writer(f)
    csv_writer.writerows([["system time", systemTime],
                          ["problemName", problemName],
                          ["N", N],
                          ["maxgen", maxgen],
                          ["D", D],
                          ["M", M]
                          ])

    while gen <= maxgen:
        print("generation for %.2f%%" % (gen / maxgen * 100))

        # recording...
        if gen % 10 == 0:
            delta = 0
            flag = 1
            if problemName == 'swimmer':
                delta = 300
                flag = -1
            csv_writer.writerows([["generation", gen]])
            for i in range(N):
                out = []
                for k in range(M):
                    out.append(flag * pop[i].obj[k] + delta)
                csv_writer.writerow(out)
            f.flush()

        pop_parent1 = tournamentSelection(pop, N)
        pop_parent2 = tournamentSelection(pop, N)
        pop_parent = pop_parent1 + pop_parent2
        pop_offspring = crossMutation(pop_parent, D, lower, upper, pc, pm, problem)
        pop_combine = pop + pop_offspring
        F, pop_combine_non, _ = nonDominatedSort(len(pop_combine), pop_combine)

        # preference
        if gen >= tau:
            # init model
            # if gen == tau:
            #     CreateModel(M, 0)
            if gen % tau == 0:
                # # collect dataset
                # number_of_inquiries = mu
                # while number_of_inquiries > 0:
                #     index = random.sample(range(0, np.shape(pop_combine_non)[0]), 2)
                #     xxx = index[0]
                #     yyy = index[1]
                #     ind1 = pop_combine_non[xxx]
                #     ind2 = pop_combine_non[yyy]
                #     if golden_function(ind1, M) > golden_function(ind2, M):
                #         winners.append(ind1.obj)
                #         losers.append(ind2.obj)
                #     else:
                #         winners.append(ind2.obj)
                #         losers.append(ind1.obj)
                #     number_of_inquiries -= 1
                # # train model
                # TrainModel(winners, losers, 0)
                # # use model
                # inputIndividual = []
                # for i in range(np.shape(pop_combine_non)[0]):
                #     inputIndividual.append(pop_combine_non[i].obj)
                # score = UseModel(inputIndividual, 0)
                for i in range(np.shape(pop_combine_non)[0]):
                    # pop_combine_non[i].cd = score[i]
                    pop_combine_non[i].cd = golden_function(pop_combine_non[i], M)
            else:
                # # use model
                # inputIndividual = []
                # for i in range(np.shape(pop_combine_non)[0]):
                #     inputIndividual.append(pop_combine_non[i].obj)
                # score = UseModel(inputIndividual, 0)
                for i in range(np.shape(pop_combine_non)[0]):
                    # pop_combine_non[i].cd = score[i]
                    # pop_combine_non[i].cd = 0.0
                    pop_combine_non[i].cd = golden_function(pop_combine_non[i], M)

        else:
            pop_combine_non = crowdingDistance(F, pop_combine_non, M)

        pop = elitism(N, pop_combine_non)
        gen = gen + 1

    end = time.time()
    print("EA runtime：%2fs" % (end - start))

    f.close()

# In[ ]:

if __name__ == "__main__":
    #oprint('1')
    main()


