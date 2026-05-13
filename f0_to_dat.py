import numpy as np 
import pyana
import sys

input = sys.argv[1]
output = sys.argv[2]

data = pyana.fzread(input)["data"]
print(data.shape)

np.savetxt(output, data[0,0], fmt="%1.5e")