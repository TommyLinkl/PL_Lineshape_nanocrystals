import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.fft import fft, ifft
from scipy.integrate import simpson

from allFunc import *

allCoupling = read_Vklq("Vklq-diabatic-new.dat", 14742, 8, 1)
g = open("Vklq-diabatic-new_trim.dat", "w")
for Ph in range(len(allCoupling)):
    g.write("%d 0 0 %.15f\n" % (Ph, allCoupling[Ph,0]))
g.close()
