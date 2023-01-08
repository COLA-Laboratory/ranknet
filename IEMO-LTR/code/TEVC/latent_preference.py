#!/usr/bin/env python
# coding: utf-8

import numpy as np
import math
import csv



def DM_surrogate_cosine(individual, M):
    with open("goldenPoints.csv","r") as file_:
        reader = csv.reader(file_)
        referencePoints = list(reader)
    referencePoints = np.array(referencePoints).astype('float32')
    referencePoint = np.reshape(referencePoints[0], (1, M))
    sum_0 = 0
    sum_1 = 0
    sum_2 = 0
    for i in range(M):
        sum_0+=individual.obj[i]*referencePoint[0][i]
        sum_1+=individual.obj[i]*individual.obj[i]
        sum_2+=referencePoint[0][i]*referencePoint[0][i]
    return sum_0/(np.sqrt(sum_1)*np.sqrt(sum_2))
    
def DM_surrogate_Tchebycheff(individual, M):
    with open("referencePoints.csv","r") as file_:
        reader = csv.reader(file_)
        referencePoints = list(reader)
    referencePoints = np.array(referencePoints).astype('float32')
    referencePoint = np.reshape(referencePoints[0], (1, M))
    sum_0 = 0
    sum_1 = 0
    sum_2 = 0
    for i in range(M):
        sum_0+=individual.obj[i]*referencePoint[0][i]
        sum_1+=individual.obj[i]*individual.obj[i]
        sum_2+=referencePoint[0][i]*referencePoint[0][i]
    return sum_0/(np.sqrt(sum_1)*np.sqrt(sum_2))


# bigger is better
def golden_function(individual, M):
    with open("weight.csv","r") as file_:
        reader = csv.reader(file_)
        weight = list(reader)
    weight = np.array(weight).astype('float32')
    weight = np.reshape(weight[0], (1, M))
    res = 0
    for i in range(M):
        tmp = abs(individual.obj[i]) / weight[0][i]
        if tmp > res:
            res = tmp

    return -res

