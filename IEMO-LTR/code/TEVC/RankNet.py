#!/usr/bin/env python
# coding: utf-8
import os
import sys

os.environ["MKL_NUM_THREADS"] = '4'
os.environ["NUMEXPR_NUM_THREADS"] = '4'
os.environ["OMP_NUM_THREADS"] = '4'
os.environ["VECLIB_MAXIMUM_THREADS"] = '4'
os.environ["OPENBLAS_NUM_THREADS"] = '4'
os.environ["TF_CPP_MIN_LOG_LEVEL"] = '2'

import numpy as np

import tensorflow as tf
from tensorflow import keras
from keras.models import Model, Sequential, load_model
from keras.layers import Dense, Input, Subtract, Activation
from keras.utils import plot_model
from keras.callbacks import TensorBoard
from keras import backend as K
from keras import initializers
import h5py

# opt para
# from sklearn.model_selection import GridSearchCV
# from keras.wrappers.scikit_learn import KerasClassifier

# create->save
def CreateModel(input_size,nthread):
    lr = 0.01  # learning rate
    layer_1 = Dense(5, activation="relu", kernel_initializer=initializers.LecunUniform(seed=None))
    layer_2 = Dense(64, activation="relu", kernel_initializer=initializers.LecunUniform(seed=None))
    layer_3 = Dense(32, activation="relu", kernel_initializer=initializers.LecunUniform(seed=None))
    s = Dense(1, name='score')

    # rel_doc -> RankNet -> rel_score
    rel_doc = Input(shape=(input_size,))
    l_1_rel = layer_1(rel_doc)
    l_2_rel = layer_2(l_1_rel)
    l_3_rel = layer_3(l_2_rel)
    rel_score = s(l_3_rel)

    # irr_doc -> RankNet -> irr_score
    irr_doc = Input(shape=(input_size,))
    l_1_irr = layer_1(irr_doc)
    l_2_irr = layer_2(l_1_irr)
    l_3_irr = layer_3(l_2_irr)
    irr_score = s(l_3_irr)

    # Subtract scores.
    diff = Subtract()([rel_score, irr_score])

    # Pass difference through sigmoid function.
    prob = Activation("sigmoid")(diff)

    # Build model.
    model = Model(inputs=[rel_doc, irr_doc], outputs=prob)
    model.compile(optimizer=keras.optimizers.Adam(learning_rate=lr),
                  loss=keras.losses.binary_crossentropy,
                  metrics=['accuracy'])

    # print('----------------CreateModel-------------------')
    # model.summary()
    #
    # # print weights
    # print('weights in layer(s) in CreateModel:')
    # for layer in model.layers:
    #     print(layer.get_weights())
    # sys.stdout.flush()

    # save the initial untrained model
    model.save('./model/initial_model_%d.h5'%nthread)
    del model


# load->train->save(update)
def TrainModel(winners_list, losers_list, nthread):
    winners = np.array(winners_list).astype('float64')
    losers = np.array(losers_list).astype('float64')
    pairsNum = winners.shape[0]
    M = winners.shape[1]
    epochs = 20
    batch_size = 5
    winners = winners.reshape(-1, M)
    losers = losers.reshape(-1, M)
    target = [1] * pairsNum
    target = np.array(target).astype('float64')
    target = target.reshape(-1, 1)

    model = load_model('./model/initial_model_%d.h5'%nthread)
    # print('----------------TrainModel-------------------')
    # model.summary()
    # print('weights in layer(s) in TrainModel(before training):')
    # for layer in model.layers:
    #     print(layer.get_weights())
    # sys.stdout.flush()

    # verbose：日志显示
    # verbose = 0 为不在标准输出流输出日志信息
    # verbose = 1 为输出进度条记录
    # verbose = 2 为每个epoch输出一行记录
    model.fit([winners, losers], target, batch_size=batch_size, epochs=epochs, verbose=0, shuffle=True)
    # print('weights in layer(s) in TrainModel(after training):')
    # for layer in model.layers:
    #     print(layer.get_weights())
    model.save('./model/model_%d.h5'%nthread)
    # print weights
    # print('weights in layer(s) after training:')
    # for layer in model.layers:
    #     print(layer.get_weights())
    del model

# load->use
def UseModel(currPop_list, nthread):
    currPop = np.array(currPop_list).astype('float64')
    currPop = currPop.reshape(-1, currPop.shape[1])

    model = load_model('./model/model_%d.h5'%nthread)
    # print('----------------UseModel-------------------')
    # model.summary()
    # print('weights in layer(s) in UseModel:')
    # for layer in model.layers:
    #     print(layer.get_weights())
    # sys.stdout.flush()

# for layer in model.layers:
    #     print('--------------------------')
    #     print(layer.get_weights())
    #     print(layer.input)
    #     print(layer.output)

    get_score = K.function([model.layers[0].input],
                           [model.get_layer('score').output])
    score = (get_score(currPop)[0]).ravel()
    score = score.tolist()
    del model
    return score