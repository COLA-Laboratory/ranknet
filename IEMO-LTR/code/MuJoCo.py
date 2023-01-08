import numpy as np
import math
import os, sys
import time
import environments
from environments.swimmer import SwimmerEnv
from environments.swimmer_v3 import SwimmerV3Env
from environments.ant import AntEnv
from environments.walker2d import Walker2dEnv
from environments.half_cheetah import  HalfCheetahEnv
from environments.hopper import HopperEnv
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
sys.path.append(base_dir)
from copy import deepcopy

def Swimmer(x):
    A = x[:-2]
    A = np.array(A).reshape(2, -1)
    b = x[16:]
    b = np.array(b).reshape(-1, 1)

    env = SwimmerEnv()
    state = env.reset()
    state = state.reshape(-1, 1)
    action = np.dot(A, state) + b

    itr = 500
    count = 0
    obj1 = 0
    obj2 = 0
    while count < itr:
        action = action.ravel()
        state, _, done, info = env.step(action)
        rewards = info["obj"]
        obj1 += rewards[0]
        obj2 += rewards[1]
        state = state.reshape(-1, 1)
        action = np.dot(A, state) + b
        count += 1

    obj = [-obj1 + 300, -obj2 + 300]
    return obj

def SwimmerV3(x):
    A = x[:-2]
    A = np.array(A).reshape(2, -1)
    b = x[16:]
    b = np.array(b).reshape(-1, 1)

    env = SwimmerV3Env()
    state = env.reset()
    state = state.reshape(-1, 1)
    action = np.dot(A, state) + b

    itr = 500
    count = 0
    obj1 = 0
    obj2 = 0
    obj3 = 0
    while count < itr:
        action = action.ravel()
        state, _, done, info = env.step(action)
        rewards = info["obj"]
        obj1 += rewards[0]
        obj2 += rewards[1]
        obj3 += rewards[2]
        state = state.reshape(-1, 1)
        action = np.dot(A, state) + b
        count += 1

    obj = [-obj1 + 300, -obj2 + 300, -obj3 + 300]
    return obj

def Hopper(x):
    print(np.shape(x))
    A = x[:-3]
    A = np.array(A).reshape(3, -1)
    print(np.shape(A))

    b = x[33:]
    b = np.array(b).reshape(-1, 1)
    print(np.shape(b))

    env = HopperEnv()
    state = env.reset()
    state = state.reshape(-1, 1)
    print(np.shape(state))
    action = np.dot(A, state) + b

    itr = 512
    count = 0
    obj1 = 0
    obj2 = 0
    while count < itr:
        action = action.ravel()
        state, _, done, info = env.step(action)
        rewards = info["obj"]
        obj1 += rewards[0]
        obj2 += rewards[1]
        state = state.reshape(-1, 1)
        action = np.dot(A, state) + b
        count += 1

    obj = [-obj1 + 6000, -obj2 + 6000]
    return obj