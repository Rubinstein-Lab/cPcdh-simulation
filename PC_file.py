# This code has been written by Nadir Boni from the Rotem Rubinster lab, Tel Aviv University, Israel.
# This code has been used for the data presented in the paper "On the Formation of Ordered Protein Assemblies in Cell-cell Interfaces" published in PNAS.
# This file run the simulations based on the inputs inserted below. Each simulation create a directory with the extracted data. 
# For analyais and visualization of the extracted data, use the 'Analysis_file.py'.



from Launcher import Initiator as init
import sys
import os
import datetime
import random


###### inputs ######

grid_len = 10  # number of squares in each side. grid size is grid_len**2. type:int
dif_trap = 100  # diffusion trap percent size to the grid. type:int, have to be between 0-100!

isoAmount = 4  # percent of proteins in the grid from all(!) isoforms. type:int, have to be between 0-100!

KAcisSelf = 8  # cis affinity to the same isoform. type:int
KAcisOther = 8  # cis affinity to other isoforms. type:int

KAtransSelf = 8  # trans affinity to the same isoform. type:int
KAtransOther = -20  # trans affinity to other isoforms. type:int

frames = 100  # number of frames. type:int
stepPerFrame = 75  # number of moves in every frame. If -1 then a constant nonuniform function will be applied. type:int
ExpNum = 1  # number of experiments - only for the data analysis function. type:int
exp_to_vis = random.randrange(1, ExpNum + 1)  # the exp number that will be visualized. If -1, all exp will be
# visualized. If 0, non will be visualized. type:int

###################


# check the input #

# check for input types #
if not isinstance(grid_len, int):
    raise IOError("grid_len is not an integer")
if not isinstance(isoAmount, int):
    raise IOError("isoAmount is not an integer")
if not isinstance(KAcisSelf, int):
    raise IOError("KAcisSelf is not an integer")
if not isinstance(KAcisOther, int):
    raise IOError("KAcisOther is not an integer")
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
    raise IOError("stepPerFrame has to be larger the -2")

# check that the diffusion trap percent is between 0 to 100 #
if dif_trap < 0 or dif_trap > 100:
    raise IOError("the diffusion trap percent is not between 0-100")

# check that the isoform amount percent is between 0 to 100 #
if isoAmount < 0 or isoAmount > 100:
    raise IOError("the isoAmount percent is not between 0-100")

######


# transform the data #

# calculate the number of proteins from each isoform based on the percentage #
input_isoamount = isoAmount
isoAmount = (isoAmount / 100) * (grid_len**2)
isoAmount = round(isoAmount)
isoAmount = isoAmount  # calculate how much from each isoform
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
name = "./Exp_"+str(input_isoamount)+"-"+str(dif_trap)+"_CIS"+str(KAcisSelf)+"_TRANS"+str(KAtransSelf)+"_Frames"+str(frames)+"_ExpNo"+str(ExpNum)
os.mkdir(name)

# create a log file with the input data in the directory #
log_file = open(name+"/info.log", "x")
log_file.write(str(frames) + "," + str(ExpNum) + "," + str(input_exp_to_vis) + "," + str(stepPerFrame) + "," + str(grid_len) + "\n" +  # for analysis later
               "# LOG FILE\n\n" +
               "date and time: " + str(datetime.datetime.now()) + "\n\n" +
               "Grid length: " + str(grid_len) + "\t" +
               "Trap area: "+str(dif_trap)+"\n" +
               "Number of proteins from each isoforms in each membrane: "+str(isoAmount)+"\n\n" +
               "Ka-Cis to self: "+str(KAcisSelf)+"\t" +
               "Ka-Cis to others: "+str(KAcisOther)+"\n" +
               "Ka-Trans to self: "+str(KAtransSelf)+"\t" +
               "Ka-Trans to others: "+str(KAtransOther)+"\n\n" +
               "Number of samples (repetitions): "+str(ExpNum)+"\t" +
               "Samples with full analytic data: "+str(exp_to_vis)+"\n" +
               "Frames: "+str(frames)+"\t" +
               "Steps between frames: "+str(stepPerFrame)+"\n\n\n")
log_file.close()

######


# send the job to initiator #

init.multi_exp(grid_len, dif_trap, isoAmount, KAcisSelf, KAcisOther, KAtransSelf, KAtransOther, frames,
               stepPerFrame, ExpNum, exp_to_vis, name)

######
