import platform
import os
import json
import copy
import math
import matplotlib.pyplot as plt

# change to directory of the script
abspath = os.path.abspath(__file__)
mydir = os.path.dirname(abspath)
os.chdir(mydir)
mydir = os.getcwd()

# create python dictionary to hold data
temp_parameters = {
    "n": 9,
    "mx": [101, 201,401,1001],
    "my": [101, 201,401,1001],
    "omega": 1,
    "mstep": 1000,
    "cas": 1,
    "fcas": 1,
    "bc": 2,
    "eos": 2,
    "dbcas": 1,
    "R": 1.0,
    "A": 1.0,
    "B": 4.0,
    "T": 0.0707475,
    "lx": [100, 200,400,1000],
    "ly": [100, 200,400,1000],
    "c": 1,
    "rho_l": [0.31, 0.34],
    "rho_g": 0.011,
    "radius": [15,20,25,30,35],
    "freq": 100,
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
if not os.path.exists("figures"):
    os.makedirs("figures")
else:
    print(f"Directory 'figures' already exists.")
os.system('g++ -o main.exe s-c-cleaned-version-reading-and-writing.cpp')
filename=[]
for i in range(len(temp_parameters["mx"])):
    parameters = copy.deepcopy(temp_parameters)
    parameters["mx"] = int(temp_parameters["mx"][i])
    parameters["my"] = int(temp_parameters["my"][i])
    parameters["lx"] = int(temp_parameters["lx"][i])
    parameters["ly"] = int(temp_parameters["ly"][i])

    for j in range(len(temp_parameters["rho_l"])):
        parameters["rho_l"] = temp_parameters["rho_l"][j]
        
        for k in range(len(temp_parameters["radius"])):
            parameters["radius"] = temp_parameters["radius"][k]
            filename.append(f"../testcase2/pressure_{parameters['lx']}_{parameters['mx']}_{parameters['rho_l']}_{parameters['radius']}_VSM.dat")
            print(f"mx = {parameters['mx']}, my = {parameters['my']}, lx = {parameters['lx']}, ly = {parameters['ly']}, rho_l = {parameters['rho_l']}, radius = {parameters['radius']}")

            # write dictionary to json file
            with open("input.json", "w") as outfile:
                json.dump(parameters, outfile)

            # windows
            # if (platform.system() == 'Windows'):
            #     os.system('main.exe')
# print(filename)
###############equation of state 
def calculate_pressure(dense,A,B,T,R):
    p=0.0
    p=dense*R*T*(1.0+B*dense/4.0+(B*dense/4.0)*(B*dense/4.0)-math.pow(B*dense/4.0,3.0))/math.pow((1.0-B*dense/4.0),3.0)-A*dense*dense
    return p   
# #######################################################
# #R-P equation of matlab code
# # create python dictionary to hold data
os.chdir(mydir+"/matlab_r_p")
rho_l=[]
for k in range(len(temp_parameters["rho_l"])):
    rho_l.append(temp_parameters["rho_l"][k])
rho_g=temp_parameters["rho_g"]
omega=temp_parameters["omega"]
lx=[]
for l in range(len(temp_parameters["lx"])):
    lx.append(temp_parameters["lx"][l])
radius=[]
for l in range(len(temp_parameters["radius"])):
    radius.append(temp_parameters["radius"][l])
A=temp_parameters["A"]
B=temp_parameters["B"]
T=temp_parameters["T"]
R=temp_parameters["R"]
p_inf=[]
for o in range(len(rho_l)):
    p_inf.append(calculate_pressure(rho_l[o],A,B,T,R))
######## variable p_v #################
def read_first_column(filename):
    P_v = []
    with open(filename, 'r') as file:
        for line in file:
            # Split the line by whitespace and take the first element
            value = float(line.strip().split()[1])
            P_v.append(value)
    return P_v
# p_v=calculate_pressure(rho_g,A,B,T,R)
# filename="/testcase2/pressure_100_101_0.31_25_VSM.dat"
# filename="/testcase2/pressure_lx_mx_rho_l_radius_VSM.dat"
p_v = [[] for _ in range(len(filename))]
# p_v=[[],[],[],[],[],[],[]]
# radius=15
k=0
for i in range(len(filename)):
    if i in [0, 5, 10, 15, 20, 25, 30, 35]:
        continue
    else:
    # print(i)
    # 15,20,25,30,35
    #15: 0,5,10,15,20,25,30,35
    #20: 1,6,11,16,21,26,31,36
    #25: 2,7,12,17,22,27,32,37
    #30: 3,8,13,18,23,28,33,38
    #35: 4,9,14,19,24,29,34,39
    # if radius==15:
    #     i=0
    # elif radius==20:
    #     i=1
    # elif radius==25:
    #     i=2
    # elif radius==30:
    #     i=3
    # elif radius==35:
    #     i=4
        # while(i<40):
        p_v[k]=read_first_column(filename[i])
        # i=i+5
        k=k+1
    # if(i==4,9,14,19,24,29):
    #     # print(i)
    # plt.plot(p_v[1], label='p_v')
######## variable p_v #################
r_inf=[]
for v in range(len(lx)):
    # r_inf.append(lx[v]/2)
    r_inf.append(lx[v]/2*0.8814)
nu_l=(1/omega-0.5)/3 #omega=1.0/(3.*alpha+0.5)
temp_parameters_rp= {
    "radius": radius,
    "P_v": p_v,
    "P_inf": p_inf,
    "S": 0.02193,
    "nu_l": nu_l,
    "rho_l": rho_l,
    "r_inf": r_inf 
}
k=0
for i in range(len(temp_parameters_rp["r_inf"])):
    parameters = copy.deepcopy(temp_parameters_rp)
    parameters["r_inf"] = temp_parameters_rp["r_inf"][i]
    for j in range(len(temp_parameters_rp["P_inf"])):
        parameters["P_inf"] = temp_parameters_rp["P_inf"][j]
        parameters["rho_l"] = temp_parameters_rp["rho_l"][j]

        for l in range(len(temp_parameters_rp["radius"])):
            # print(l)
            if l==0:
                continue
            else:
                parameters["radius"] = temp_parameters_rp["radius"][l]
                parameters["P_v"] = temp_parameters_rp["P_v"][k]
                print(f"P_inf = {parameters['P_inf']}, rho_l = {parameters['rho_l']}, r_inf = {parameters['r_inf']},radius = {parameters['radius']}")
    # write dictionary to json file
                with open(f"input{k+1}.json", "w") as outfile:
                    json.dump(parameters, outfile)
                k=k+1
# start matlab terminal and run test case
# os.system('matlab -nodesktop -nosplash -r "dr_post_measure_difference; exit;"')


