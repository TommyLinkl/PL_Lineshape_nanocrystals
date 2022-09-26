import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

allCoupling = read_Vklq("Vklq_diabatic.dat", 7182, 20, 1)

for exc_index in range(len(allCoupling[0])): 
    g = open("Vklq_trim_Exc%d.dat" % exc_index, "w")
    for Ph in range(len(allCoupling)):
        g.write("%d %d %d %.15f\n" % (Ph, exc_index, exc_index, allCoupling[Ph,exc_index]))
    g.close()
