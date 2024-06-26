# This code has been written by Nadir Boni from the Rotem Rubinster lab, Tel Aviv University, Israel.
# This code has been used for the data presented in the paper "On the Formation of Ordered Protein Assemblies in Cell-cell Interfaces" published in PNAS.
# This file is being accesed by the 'PC_file.py' and responsible for the initation of all the exprement as instructed by the 'PC_file.py' code.



import os
import time
import random
from Launcher import Functions as fn
import zipfile


# main function that conduct the experiments one by one and save the data to the directory. This function only create
# the data files, analysis need to be conducted later with the envelope_analysis.py file #

def multi_exp(grid_len, dif_trap, isoAmount, KAcisSelf, KAcisOther, KAtransSelf, KAtransOther, frames,
              stepPerFrame, ExpNum, exp_to_vis, name):
    t1 = time.time()  # start time for running time calculations

    nonuniform_flag = False
    if stepPerFrame == -1:  # for nonuniform steps per frame
        nonuniform_flag = True

    for exp in range(1, ExpNum + 1):
        print(exp)
        # create the initial grid with cis-dimers, using the Functions.py file
        grid, all_protein, all_complex = fn.init_grid(grid_len, isoAmount)
        trap_locations = fn.dif_trap_loc(grid_len, dif_trap)

        # creating a sub-directory for the experiment with two data files (one with general stats in each frame and
        # second with the length of all complexes in each frame)
        os.mkdir(name + "/Exp" + str(exp))
        filename1 = name + "/Exp" + str(exp) + "/Stats_Exp" + str(exp)
        file1 = open(filename1, 'w')
        file1.write("# General stats about each frame in the exp. Each line represent a frame #\n" +
                    "# num of cis, num of trans, avg complex length, max complex length\n")  # add header to the file
        file1.close()

        filename2 = name + "/Exp" + str(exp) + "/Full_Complex_Length_Exp" + str(exp)
        file2 = open(filename2, 'w')
        file2.close()

        # create sub-directory for raw data for later analysis
        if exp in exp_to_vis:
            os.mkdir(name+"/Exp"+str(exp)+"/Raw_data")

        for frame in range(frames):

            # check that proteins haven't overridden between frames
            if sum(sum(sum(grid))) != isoAmount * 5:
                print("wrong here")
                fn.save_grid(all_protein, all_complex, exp, frame, name)
                raise SystemError("protein lost")

            # check that all protein have the same direction as their complex
            for p in all_protein:
                if p.direction not in p.complex.direction:
                    print(grid)
                    print("wrong there")
                    print(p.complex.proteins)
                    raise SystemError("protein in wrong direction")

            # if the current running experiment is for visualization full grid data is being saved
            if exp in exp_to_vis:
                fn.save_grid(all_protein, all_complex, exp, frame, name)

            # update the two data files after each frame
            fn.analyse_grid(filename1, filename2, all_protein, all_complex)

            if nonuniform_flag:
                stepPerFrame = fn.StepNonUni(frame)

            for step in range(stepPerFrame):
                what_move = random.choice(['move', 'turn', 'create_cis', 'break_cis', 'create_trans', 'break_trans'])

                if what_move == 'move':
                    grid = fn.move(grid, all_protein, trap_locations)
                elif what_move == 'turn':
                    grid = fn.turn(grid, all_protein, trap_locations)
                elif what_move == 'create_cis':
                    grid, all_complex, all_protein = fn.create_cis(grid, all_protein, all_complex, KAcisSelf,
                                                                   KAcisOther)
                elif what_move == 'create_trans':
                    grid, all_complex = fn.create_trans(grid, all_protein, all_complex, KAtransSelf, KAtransOther,
                                                        trap_locations)
                elif what_move == 'break_cis':
                    grid, all_complex, all_protein = fn.break_cis(grid, all_protein, all_complex, KAcisSelf, KAcisOther,
                                                                  frame, step)
                elif what_move == 'break_trans':
                    grid, all_complex = fn.break_trans(grid, all_protein, all_complex, KAtransSelf, KAtransOther, frame,
                                                       step)

        # update the two data files in the last frame
        if exp in exp_to_vis:
            fn.save_grid(all_protein, all_complex, exp, frame+1, name)

            # zipping the "Raw_data" directory to save memory
            zip_file = zipfile.ZipFile(name + "/Exp" + str(exp) + "/Raw_data.zip", 'w')
            for dirname, subdirs, files in os.walk(name + "/Exp" + str(exp) + "/Raw_data"):
                for file in files:
                    filepath = os.path.join(dirname, file)
                    zip_file.write(filepath, os.path.basename(filepath))  # zip each file
                    os.remove(filepath)  # remove the original file
            zip_file.close()
            os.rmdir(name + "/Exp" + str(exp) + "/Raw_data")  # remove the original directory

        fn.analyse_grid(filename1, filename2, all_protein, all_complex)

    t2 = time.time() # finish time for running time calculations

    # updating the log file
    log_file = open(name + "/info.log", "a")
    log_file.write("\nStarting to run experiments...\n" +
                   "\tStar time: " + str(t1) + "\n" +
                   "\tConducted experiments: " + str(ExpNum) + "\n" +
                   "\tTrap locations in the grid: " + str(trap_locations) + "\n" +
                   "\tFinish time: " + str(t2) + "\n" +
                   "\tRun time: " + str((t2 - t1)/60) + " minutes\n" +
                   "Finished conducting experiments!\n")
    log_file.close()
    print(trap_locations)

    print("run time:", (t2 - t1) / 60)
    return None


