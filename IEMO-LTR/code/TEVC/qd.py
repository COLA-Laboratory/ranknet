import math
import os, sys
import time
import environments
from environments.swimmer import SwimmerEnv
from environments.ant import AntEnv
from environments.walker2d import Walker2dEnv
from environments.half_cheetah import  HalfCheetahEnv
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
sys.path.append(base_dir)
import numpy as np
from copy import deepcopy

def main():
    env = SwimmerEnv()
    # envtemp = deepcopy(env)


    count=0
    a1 = 0.5
    a2 = 0.5
    for i in range(100):

        #env1= deepcopy(env)
        state = env.reset_model()
        count = 0
        while count <1:
            action = np.array([a1,a2])
            state,_,done,info = env.step(action)
            print(info)
            count += 1


if __name__ == "__main__":
    #oprint('1')
    main()







