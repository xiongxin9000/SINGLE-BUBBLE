import platform
import os
import json
import copy
import math
import matplotlib.pyplot as plt
import numpy as np
rho_liquid=[]
rho_gas=[]
# change to directory of the script
abspath = os.path.abspath(__file__)
mydir = os.path.dirname(abspath)
os.chdir(mydir)
mydir = os.getcwd()
# iR_values = []
p_values = []
T_reduced = [0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.970339,0.988861,1.0]
T_ref=0.09433
T_lbm=[]
rho_liquid_temp=[]
rho_gas_temp=[]
#digitize the rho values from literature and set the initial rho_l and rho_g for LBM to simulate 
#after it reaches the stable state, take this as an equilibrium value for rho_l and rho_g in LBM
for T in T_reduced:
    T_lbm.append(T*T_ref)
    #RESERVE 6 DIGITS OF PRECISION
    T_lbm = [round(x,6) for x in T_lbm]
print(T_lbm)
# create python dictionary to hold data
temp_parameters = {
    "n": 9,
    "mx": 21,
    "my": 201,
    "omega": 1,
    "mstep": 10000,
    "cas": 10,
    "fcas": 1,
    "bc": 1,
    "eos": 2,
    "dbcas": 1,
    "R": 1.0,
    "A": 1.0,
    "B": 4.0,
    "T": T_lbm,#different T
    "lx": 20,
    "ly": 200,
    "c": 1,
    "rho_l": [0.43,0.41,0.38,0.36,0.33,0.31,0.28,0.25,0.21,0.19,0.16,0.13],#different rho_l
    "rho_g": [0.0002,0.00046,0.0021,0.0055,0.011,0.019,0.029,0.044,0.066,0.08,0.1,0.13],#different rho_g
    "radius": 20,
    "freq": 1000,
    "rho0": 1.0,
    "g": -1,
    "if_th": 5.0,
    "path": "./testcase2/"
}
with open("input_matlab.json", "w") as outfile:
    json.dump(temp_parameters, outfile)
#create a folder for output
if not os.path.exists(temp_parameters['path']):
    os.makedirs(temp_parameters['path'])
else:
    print(f"Directory '{temp_parameters['path']}' already exists.")
#create a folder for figures
if not os.path.exists("figures2"):
    os.makedirs("figures2")
else:
    print(f"Directory 'figures2' already exists.")
os.system('g++ s-c-cleaned-version-reading-and-writing.cpp -o s-c-cleaned-version-reading-and-writing')

for i in range(len(temp_parameters["T"])):
    parameters = copy.deepcopy(temp_parameters)
    parameters["T"] = temp_parameters["T"][i]
    parameters["rho_l"] = temp_parameters["rho_l"][i]
    parameters["rho_g"] = temp_parameters["rho_g"][i]
    print(f"mx = {parameters['mx']}, my = {parameters['my']}, \
          rho_l = {parameters['rho_l']}, rho_g = {parameters['rho_g']},\
          T = {parameters['T']}")

    # write dictionary to json file
    with open("input.json", "w") as outfile:
        json.dump(parameters, outfile)

    # windows
    if (platform.system() == 'Windows'):
        os.system('s-c-cleaned-version-reading-and-writing.exe')
###############maxwell-area-construction#################################
###################python-version#####################
#create an array to store rho_l and rho_g with different temperature
    data = np.loadtxt(f'testcase2/{str(parameters["rho_l"])}.dat')
    rho1=data[-1,0]
    rho2=data[-1,1]
    rho_liquid_temp.append(rho1)
    rho_gas_temp.append(rho2)
rho_liquid=np.append([],rho_liquid_temp)
rho_gas=np.append([],rho_gas_temp)
# plot the data
plt.plot(rho_liquid,T_reduced)
plt.plot(rho_gas,T_reduced)
plt.xlabel('density')
plt.ylabel('T')
plt.title('maxwell-area-construction')
# plt.show()
plt.savefig(f'figures2/maxwell-area-construction.png')
################matlab-version################################
# # start matlab terminal and run test case
# os.system('matlab -nodesktop -nosplash -r "maxwell; exit;"')
