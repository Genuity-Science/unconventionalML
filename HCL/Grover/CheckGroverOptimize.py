# native libraries
from GroverOptimize import *
import matplotlib.pyplot as plt
import pandas as pd
import feather

def checkSimpleOptimize():
    from GroverOptimize import *
    groverminimize(opt_oracle_lessthan, 5)
    x = [random.uniform(-1, 1) for i in range(0, 2 ** 5)]
    f = lambda i: x[i]


