import platform
import os
import json
import copy
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
# change to directory of the script
abspath = os.path.abspath(__file__)
mydir = os.path.dirname(abspath)
os.chdir(mydir)
mydir = os.getcwd()
iR_values = []
p_values = []
# create python dictionary to hold data
temp_parameters = {
    "n": 9,
    "mx": 101,
    "my": 101,
    "omega": 1,
    "mstep": 10000,
    "cas": 1,
    "fcas": 1,
    "bc": 1,
    "eos": 2,
    "dbcas": 1,
    "R": 1.0,
    "A": 1.0,
    "B": 4.0,
    "T": 0.0707475,
    "lx": 100,
    "ly": 100,
    "c": 1,
    "rho_l": 0.33,
    "rho_g": 0.011,
    "radius": [20,25,30],
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
os.system('g++ -o main.exe s-c-cleaned-version-reading-and-writing.cpp')

for i in range(len(temp_parameters["radius"])):
    parameters = copy.deepcopy(temp_parameters)
    parameters["radius"] = int(temp_parameters["radius"][i])

    print(f"mx = {parameters['mx']}, my = {parameters['my']}, lx = {parameters['lx']}, ly = {parameters['ly']}, rho_l = {parameters['rho_l']},radius = {parameters['radius']}")

    # write dictionary to json file
    with open("input.json", "w") as outfile:
        json.dump(parameters, outfile)

    # windows
    # if (platform.system() == 'Windows'):
    #     os.system('main.exe')
###############verify laplace law#################################
###################python-version#####################
#     fcas_str=""
#     if(parameters["fcas"]==1):
#         fcas_str="VSM"
#     elif(parameters["fcas"]==2):
#         fcas_str="EDM"
#     elif(parameters["fcas"]==3):
#         fcas_str="Guo"
#     data = np.loadtxt(f'testcase2/pressure_{str(parameters["lx"])}_{str(parameters["mx"])}_{str(parameters["rho_l"])}_{str(parameters["radius"])}_{fcas_str}.dat')

#     inverseR = data[:, 1]
#     pressure = data[:, 0]
#     plt.plot(inverseR)
#     plt.xlabel('t')
#     plt.ylabel('1/R')
#     plt.title('REVOLUTION OF BUBBLE')
#     # plt.show()
#     plt.savefig(f'figures2/{str(parameters["radius"])}.png')
#     plt.clf()  # Clear the current figure
# #when delta is less than 0.0001 then the system is stable 
# #take that average of last 100 steps of R and pressure to a two column array 
#     last_100_steps_inverseR = inverseR[-100:]
#     last_100_steps_pressure = pressure[-100:]
#     average_inverseR = np.mean(last_100_steps_inverseR)
#     average_pressure = np.mean(last_100_steps_pressure)
#     iR_values.append(average_inverseR)
#     p_values.append(average_pressure)
# iR = np.append([], iR_values)
# p = np.append([], p_values)
# slope, intercept, r_value, p_value, std_err = linregress(iR, p)
# print(slope)
# plt.figure(figsize=(10, 6)),
# #plot that array to verify that it is follow the laplace law.
# plt.plot(iR,p)
# plt.xlabel('1/R')
# plt.ylabel('pressure difference')
# plt.title('laplace law')
# # plt.show()
# plt.savefig('figures2/laplacelaw.png')
# # create python dictionary to hold data
################matlab-version################################
# # start matlab terminal and run test case
# os.system('matlab -nodesktop -nosplash -r "laplacelaw; exit;"')
