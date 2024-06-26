import zipfile
import os

from Launcher import Analysis as an

###### inputs #####

file_path = "./Exp"  # path to the general folder, depending on the analysis type. type:string

# the path needs to be a specific condition folder (one log file) #
Basic_graph_one_exp = True  # If True will create basic graphs and movie for exps with full data. type:bool
visualize_by = 'm'  # 't' for visualize by complex ; 'i' by isoform type ; 't' by trans interactions ; 'm' by membrane;
# 'cl' by cluster


# the path needs to be to a folder containing all 0-8 combination (81 log files) #
# for correct representation make sure all steps per frame on those exps are the same! #
Create_data_for_R = False  # If True will create heatmaps and additional data for R plotting. type:bool
conc = 4  # need to supply the constants of all 81 exps (easy to read from folder names)
trap = 0.05
expnum = 5
frames = 100
same_isoform_proportion = '3~5'  # number of same isoforms out of the total isoform in each membrane (format: 'same~all'). type:string
Annot_HeatMap = True  # add numbers values to the heatmap cells?
rng = range(6, 11, 2)  # the Ka range of experiments we want to visualize as heat maps. type: range.

###################


if Basic_graph_one_exp:
    # get the number of frames, samples, sample to vis and grid length from the log file #
    f = open(str(file_path) + "/info.log", 'r')
    line = f.readline()  # reading the first line in the file that holds the information
    line = line.strip()
    line_lst = line.split(sep=",")
    frames = int(line_lst[0])
    ExpNum = int(line_lst[1])
    exp_to_vis = int(line_lst[2])
    grid_len = int(line_lst[4])
    f.close()

    # basic analysis and movie making #
    an.graph_stat(ExpNum, frames, file_path)  # basic analysis function

    if exp_to_vis == -1:  # visualization function
        for exp in range(1, ExpNum + 1):
            # unzipping the "Raw_data" directory if haven't done before
            if not os.path.isdir(str(file_path) + "/Exp" + str(exp) + "/Raw_data"):
                os.mkdir(str(file_path) + "/Exp" + str(exp) + "/Raw_data")
                with zipfile.ZipFile(str(file_path) + "/Exp" + str(exp) + "/Raw_data.zip", 'r') as zip_ref:
                    zip_ref.extractall(str(file_path) + "/Exp" + str(exp) + "/Raw_data")
            print("Unzipped exp " + str(exp))
            # analysing the data as instructed by the inputs
            if visualize_by == 'i':
                an.vis_by_interaction(exp, frames, grid_len, file_path)
            elif visualize_by == 'c':
                an.vis_by_complex(exp, frames, grid_len, file_path)
            elif visualize_by == 't':
                an.vis_by_trans(exp, frames, grid_len, file_path)
            elif visualize_by == 'm':
                an.vis_by_membrane(exp, frames, grid_len, file_path)
            elif visualize_by == 'cl':
                an.vis_by_cluster(exp, frames, grid_len, file_path)
            else:
                raise IOError("visualize_by input is unknown")
            # creating the movie
            an.mve_from_img(exp, frames, file_path)
    elif exp_to_vis == 0:
        pass
    else:
        # unzipping the "Raw_data" directory if haven't done before
        if not os.path.isdir(str(file_path) + "/Exp" + str(exp_to_vis) + "/Raw_data"):
            os.mkdir(str(file_path) + "/Exp" + str(exp_to_vis) + "/Raw_data")
            with zipfile.ZipFile(str(file_path) + "/Exp" + str(exp_to_vis) + "/Raw_data.zip", 'r') as zip_ref:
                zip_ref.extractall(str(file_path) + "/Exp" + str(exp_to_vis) + "/Raw_data")
        print("Unzipped exp " + str(exp_to_vis))
        # analysing the data
        if visualize_by == 'i':
            an.vis_by_interaction(exp_to_vis, frames, grid_len, file_path)
        elif visualize_by == 'c':
            an.vis_by_complex(exp_to_vis, frames, grid_len, file_path)
        elif visualize_by == 't':
            an.vis_by_trans(exp_to_vis, frames, grid_len, file_path)
        elif visualize_by == 'm':
            an.vis_by_membrane(exp_to_vis, frames, grid_len, file_path)
        elif visualize_by == 'cl':
            an.vis_by_cluster(exp_to_vis, frames, grid_len, file_path)
        else:
            raise IOError("visualize_by input is unknown")
        # creating the movie
        an.mve_from_img(exp_to_vis, frames, file_path)

    print("Done basic data for one exp")

if Create_data_for_R:
    # check that all steps per frame is the same #
    old_steps = None  # to keep a reference value for later comparison
    for i in rng:  # cis = i, trans = j
        for j in rng:
            f = open(
                file_path + "/Exp_" + str(conc) + "-" + str(trap) + "_CIS" + str(i) + "_TRANS" + str(j) + "_Frames" +
                str(frames) + "_ExpNo" + str(expnum) + "/info.log", "r")
            line = f.readline()
            line = line.strip()
            line_lst = line.split(sep=",")
            new_steps = int(line_lst[3])  # extracting the 'steps per frame' data from the first line in the log file
            if old_steps is None:
                old_steps = int(line_lst[3])  # updating 'old_steps' for the first time
            else:
                if old_steps != new_steps:
                    raise IOError("The number of steps in all the exps is different")
                else:
                    continue
    print("Finished 'steps per frame' checking process")

    # unzipping "Raw_data" directory if it hasn't unzipped before
    for cis in rng:  # cis = i, trans = j
        for trans in rng:
            file_path_2 = str(file_path) + "/Exp_" + str(conc) + "-" + str(trap) + "_CIS" + str(cis) + "_TRANS" + \
                          str(trans) + "_Frames" + str(frames) + "_ExpNo" + str(expnum)
            f = open(str(file_path_2) + "/info.log", 'r')
            line = f.readline()
            line = line.strip()
            line_lst = line.split(sep=",")
            exp_to_vis = int(line_lst[2])
            f.close()
            if exp_to_vis == 0:
                continue
            if not os.path.isdir(str(file_path_2) + "/Exp" + str(exp_to_vis) + "/Raw_data"):
                os.mkdir(str(file_path_2) + "/Exp" + str(exp_to_vis) + "/Raw_data")
                with zipfile.ZipFile(str(file_path_2) + "/Exp" + str(exp_to_vis) + "/Raw_data.zip", 'r') as zip_ref:
                    zip_ref.extractall(str(file_path_2) + "/Exp" + str(exp_to_vis) + "/Raw_data")
    print("Finished unzipping process")

    # analysing the data
    an.Ka_heatmap(conc, trap, expnum, frames, same_isoform_proportion, file_path, Annot_HeatMap, rng)
    an.Ka_lineplot(0, 8, file_path)  # the constant cis is set to 0 and the constant trans is set to 8
    # an.Histogram_db_creator(conc, trap, expnum, frames, file_path)
