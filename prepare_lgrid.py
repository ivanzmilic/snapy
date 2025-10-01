import numpy as np
import sys

# Here you can make some wavelength grid by hand: 

# Ca K 
ll0 = 3933.6
llrange = 10.0
Deltall = 0.01
NL = int(llrange / Deltall) + 1

ll = np.linspace(ll0 - llrange*0.5, ll0 + llrange*0.5, NL)

lgrid = ll

# Mg I  b 

ll0 = 5172.0
llrange = 5.0
Deltall = 0.01
NL = int(llrange / Deltall) + 1

ll = np.linspace(ll0 - llrange*0.5, ll0 + llrange*0.5, NL)

lgrid = np.append(lgrid, ll)

# Fe 5250 

ll0 = 5250.2
llrange = 2.0
Deltall = 0.01
NL = int(llrange / Deltall) + 1

ll = np.linspace(ll0 - llrange*0.5, ll0 + llrange*0.5, NL)

lgrid = np.append(lgrid, ll)
# Fe Na I D

ll0 = 5893.0
llrange = 10.0
Deltall = 0.01
NL = int(llrange / Deltall) + 1

ll = np.linspace(ll0 - llrange*0.5, ll0 + llrange*0.5, NL)

lgrid = np.append(lgrid, ll)

# Ca II 8542

ll0 = 8542.0
llrange = 20.0
Deltall = 0.02
NL = int(llrange / Deltall) + 1

ll = np.linspace(ll0 - llrange*0.5, ll0 + llrange*0.5, NL)

lgrid = np.append(lgrid, ll)

NLtot = len(lgrid)

print("info::total number of points would be: ", NLtot)

np.savetxt("acoustic_lgrid.dat", lgrid.T, header = str(NLtot), comments = '')
