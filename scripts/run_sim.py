import numpy as np
import matplotlib.pyplot as plt
from src.worldModels import *
import scipy.io
import os.path
from src.utils import *
from src.run_simulations import *
from src.decoding import *
from src.agents import *
from sklearn import svm
import pickle
from tqdm import tqdm
import ipywidgets as widgets
from ipywidgets import interact
from IPython.display import display
from sklearn.linear_model import LinearRegression
from src.expfitting import *


rlow = 0
metrics = get_IB_metrics(0.4, 0.9, rlow=rlow)
print(metrics)