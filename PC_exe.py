# This code has been written by Nadir Boni from the Rotem Rubinstein lab, Tel-Aviv University, Israel.
# This code has been submitted in the paper "clustered protocadherin cis-interactions are required
# for combinatorial cell–cell recognition underlying neuronal self-avoidance" (PNAS 2024).


from Launcher import Initiator as init
import pandas as pd
import numpy as np
import os
import datetime
import random


###### inputs ######

grid_len = 50  # number of squares in each side. grid size is grid_len**2. type:int
dif_trap = 5  # diffusion trap percent size to the grid. type:int, have to be between 0-100!

num_SameIso = 1  # number of the same isoforms between the membranes. type:int
num_DifIso = 0  # number of the different isoforms between the membranes. type:int
isoAmount = 4  # percent of proteins in the grid from all(!) isoforms. type:int, have to be between 0-100!

KAcis = 6  # cis affinity. can be type:int and the cis affinity will be the same for all isoforms, and can be type:str
# and contain a path to a csv file with an affinity table

KAtransSelf = 6  # trans affinity to the same isoform. type:int
KAtransOther = -20  # trans affinity to other isoforms. type:int

frames = 100  # number of frames. type:int
stepPerFrame = 300000  # number of moves in every frame. If -1 then a nonuniform function will be applied. type:int
ExpNum = 1  # number of experiments - only for the data analysis function. type:int
exp_to_vis = random.randrange(1, ExpNum + 1)  # the exp number that will be visualized. If -1, all exp will be
# visualized. If 0, non will be visualized. type:int

###################


# check the input #

# check for input types #
if not isinstance(grid_len, int):
    raise IOError("grid_len is not an integer")
if not isinstance(num_SameIso, int):
    raise IOError("num_SameIso is not an integer")
if not isinstance(num_DifIso, int):
    raise IOError("num_DifIso is not an integer")
if not isinstance(isoAmount, int):
    raise IOError("isoAmount is not an integer")
if not isinstance(KAcis, int):
    if not isinstance(KAcis, str):
        raise IOError("KAcis is not an integer or sting")
if not isinstance(KAtransSelf, int):
    raise IOError("KAtransSelf is not an integer")
if not isinstance(KAtransOther, int):
    raise IOError("KAtransOther is not an integer")
if not isinstance(ExpNum, int):
    raise IOError("ExpNum has to be an int")
if not isinstance(exp_to_vis, int):
    raise IOError("exp_to_vis has to be an integer")
if not isinstance(dif_trap, float):
    if not isinstance(dif_trap, int):
        raise IOError("dif_trap is not a float/int")
if not isinstance(frames, int):
    raise IOError("frames is not an integer")
if not isinstance(stepPerFrame, int):
    raise IOError("stepPerFrame is not an integer")

# check that the exp we want to make a movie out of is within the range of exps we are conducting #
if exp_to_vis > ExpNum or exp_to_vis < -1:
    raise IOError("exp_to_vis has to be smaller or equal to ExpNum")

# check that the steps per frame is within the range of -1 to inf #
if stepPerFrame < -1:
    raise IOError("stepPerFrame has to be >= -1")

# check that the diffusion trap percent is between 0 to 100 #
if dif_trap < 0 or dif_trap > 100:
    raise IOError("the diffusion trap percent is not between 0-100")

# check that the isoform amount percent is between 0 to 100 #
if isoAmount < 0 or isoAmount > 100:
    raise IOError("the isoAmount percent is not between 0-100")

######


# transform the data #

# calculate the total isoforms amount in each membrane #
num_TotIso = num_SameIso + num_DifIso

# calculate the total isoforms amount in both(!) membranes #
TotIso = num_SameIso + (num_DifIso*2)

# create an array of cis affinities from file path or number as inputted #
if isinstance(KAcis, str):
    Cis_Mrx = pd.read_excel(KAcis, index_col=0)
else:
    Cis_Mrx = pd.DataFrame(np.zeros((TotIso, TotIso), int))
    Cis_Mrx = Cis_Mrx.replace(0, KAcis)

# check that the number of cells in the affinity table is sufficient #
if Cis_Mrx.shape[0] < TotIso or Cis_Mrx.shape[1] < TotIso:
    raise IOError("The affinity matrix dose not cover the amount of isoforms inputted")

print(Cis_Mrx)

# calculate the number of proteins from each isoform based on the percentage #
input_isoamount = isoAmount
isoAmount = (isoAmount / 100) * (grid_len**2)
isoAmount = round(isoAmount)
isoAmount = isoAmount //num_TotIso  # calculate how much from each isoform
if isoAmount % 2 != 0:  # check that the number of isoforms is even, since the initial grid is cis-dimers
    isoAmount = isoAmount - 1

# turn trap size from percentage to decimal #
dif_trap = dif_trap / 100

# create a list of exp to visualized for the -1, 0 options (visualize all or non) #
input_exp_to_vis = exp_to_vis
if exp_to_vis == -1:
    exp_to_vis = list(range(1, ExpNum + 1))
elif exp_to_vis == 0:
    exp_to_vis = list()
else:
    exp_to_vis = [exp_to_vis]

# create a directory for the job #
if isinstance(KAcis, int):
    name = "./Exp_" + str(input_isoamount) + "-" + str(dif_trap) + "_CIS" + str(KAcis) + "_TRANS" + str(KAtransSelf) + \
           "_Frames" + str(frames) + "_ExpNo" + str(ExpNum) + "_Identity" + \
           str(int(num_SameIso / (num_SameIso + num_DifIso)))
else:
    name = "./Exp_" + str(input_isoamount) + "-" + str(dif_trap) + "_CISmatrix_TRANS" + str(KAtransSelf) + \
           "_Frames" + str(frames) + "_ExpNo" + str(ExpNum) + "_Identity" + \
           str(int(num_SameIso / (num_SameIso + num_DifIso)))

os.mkdir(name)

# create a log file with the input data in the directory #
log_file = open(name+"/info.log", "x")
log_file.write(str(frames) + "," + str(ExpNum) + "," + str(input_exp_to_vis) + "," + str(stepPerFrame) + "," +
               str(grid_len) + "\n" +  # first line saves ordinal data for analysis functions later
               "# LOG FILE\n\n" +
               "date and time: " + str(datetime.datetime.now()) + "\n\n" +
               "Grid length: " + str(grid_len) + "\t" +
               "Trap area: " + str(dif_trap)+"\n" +
               "Number of same isoforms: " + str(num_SameIso) + "\t" +
               "Number of different isoforms: " + str(num_DifIso) + "\t" +
               "Total number of isoforms in each membrane: " + str(num_TotIso) + "\n" +
               "Number of proteins from each isoforms in each membrane: " + str(isoAmount)+"\n\n" +
               "Ka-Cis:\n" + str(Cis_Mrx)+"\n\n" +
               "Ka-Trans to self: " + str(KAtransSelf)+"\t" +
               "Ka-Trans to others: " + str(KAtransOther)+"\n\n" +
               "Number of samples (repetitions): " + str(ExpNum)+"\t" +
               "Samples with full analytic data: " + str(exp_to_vis)+"\n" +
               "Frames: " + str(frames) + "\t" +
               "Steps between frames: " + str(stepPerFrame)+"\n\n\n")
log_file.close()

######


# send the job to initiator #

init.multi_exp(grid_len, dif_trap, num_SameIso, num_DifIso, num_TotIso, isoAmount, Cis_Mrx, KAtransSelf, KAtransOther,
               frames, stepPerFrame, ExpNum, exp_to_vis, name)

######
