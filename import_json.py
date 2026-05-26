import json
import numpy as np 

# Define the file name
filename = 'config_lw_falc.json'

# Open the file and load its contents into a Python dictionary
with open(filename, 'r') as file:
    config_data = json.load(file)

# Print the top-level keys to verify it loaded correctly
print("Successfully loaded JSON. Top-level keys:")
print(config_data.keys())

# Check the atmosphere:

if 'atmosphere' in config_data:
    print("\nAtmosphere settings:")
    for key, value in config_data['atmosphere'].items():
        print(f"{key}: {value}")
else:    print("\nNo 'atmosphere' key found in the JSON data.")

# Then extract the quantities and pack into a snapi - style atmosphere:

atmos_snapi = np.zeros([12, len(config_data['atmosphere']['zgrid'])])
ND = len(config_data['atmosphere']['zgrid'])
print (atmos_snapi.shape)

atmos_snapi[0,:] = np.linspace(-8,1,ND) # log tau, just a placeholder for now
atmos_snapi[1,:] = np.asarray(config_data['atmosphere']['zgrid'])* 1E5 # convert from Km to cm
atmos_snapi[2,:] = config_data['atmosphere']['temp']
atmos_snapi[3,:] = config_data['atmosphere']['pg']
atmos_snapi[4,:] = config_data['atmosphere']['pel']
atmos_snapi[8,:] = np.ones(ND) * 3E5 # microturbulence, just a placeholder for now

np.savetxt("falc_82.dat", atmos_snapi.T, header=str(ND)+' falc82', fmt="%1.5e")