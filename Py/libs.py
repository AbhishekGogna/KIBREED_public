# =============================================================================
# MODULES AND ENVIRONMENT VARIABLES
# =============================================================================

# Python built-in
import os
import sys
import time
import logging
import json
import pickle
import re
import math
import importlib
import random
from datetime import datetime
from pprint import pformat
from os.path import exists

# Data wrangling
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_rows', 10)
import pyarrow.feather as ft
import numpy as np

# Predictions
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import pearsonr
import tensorflow as tf
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras import Model
from tensorflow.keras.layers import Dense, Dropout, Conv1D, add, Flatten, LSTM, AveragePooling1D, concatenate, MaxPool1D, BatchNormalization, Reshape, Activation
from tensorflow.keras.callbacks import EarlyStopping, Callback, ModelCheckpoint, ReduceLROnPlateau, TensorBoard
from tensorflow.keras.regularizers import l1, l2, L1L2
from tensorflow.keras.optimizers import Adam
from tensorflow.keras import backend as K
from tensorflow.keras import Input
import keras_tuner as kt
