import glob

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moviepy.editor import *
import os
from PIL import Image


# create a visualization of the grid for each of the frames by interaction type (later those images will turn into a
# movie) #
def vis_by_interaction(exp, frames, grid_len, path):
    square_pix_size = 10  # one square in the grid will be square_pix_size x square_pix_size pixels
    size = (square_pix_size * grid_len) + 1  # calculate the image size
    color_dict = {'o1': (241, 169, 167), 't1': (225, 66, 61),  'c1': (123, 23, 20),
                  'o2': (177, 213, 251), 't2': (11, 109, 213), 'c2': (8, 79, 155),
                  'o3': (250, 247, 158), 't3': (244, 237, 42), 'c3': (194, 188, 10),
                  'o4': (160, 228, 182), 't4': (80, 206, 120), 'c4': (31, 111, 57)}  # to translate the colour coded words to RGB colours
    for frame in range(frames + 1):  # create an image for each frame based on the grid data in Raw_data folder -
        # only to exp that was marked for visualization.
        img = Image.new('RGB', (size, size), color=(255, 255, 255))
        pix = img.load()
        f = open(path + "/Exp" + str(exp) + "/Raw_data/grid_frame" + str(frame), 'r')
        f.readline()
        f.readline()  # to skip the first 2 header lines
        for line in f:
            line = line.strip()
            line_lst = line.split(sep=",")
            color = color_dict[line_lst[2]]
            x = int(line_lst[0]) * square_pix_size
            y = int(line_lst[1]) * square_pix_size
            for i in range(square_pix_size):
                for j in range(square_pix_size):
                    pix[x + i, y + j] = color
        for i in range(0, size, square_pix_size):  # to create white spaces between each square in the grid
            for j in range(size):
                pix[i, j] = (255, 255, 255)
                pix[j, i] = (255, 255, 255)
        f.close()

        # correct the image format so it will be ready for mp4 format
        img = img.resize((1000, 1000))
        img = img.convert("RGBA")

        # save the image
        img.save(path + "/Exp" + str(exp) + "/Raw_data/frame" + str(frame) + ".png")


# create a visualization of the grid for each of the frames by complexID (later those images will turn into a
# movie) #
def vis_by_complex(exp, frames, grid_len, path):
    square_pix_size = 10  # one square in the grid will be square_pix_size x square_pix_size pixels
    size = (square_pix_size * grid_len) + 1  # calculate the image size
    color_lst = [(126, 160, 129), (109, 174, 219), (225, 187, 201), (203, 1, 21), (255, 123, 71), (219, 212, 107),
                 (70, 70, 53), (115, 89, 102), (7, 157, 109), (14, 73, 88)]  # to translate the colour coded ID
    # complex to RGB colours

    for frame in range(frames + 1):  # create an image for each frame based on the grid data in Raw_data folder -
        # only to exp that was marked for visualization.
        img = Image.new('RGB', (size, size), color=(255, 255, 255))
        pix = img.load()
        f = open(path + "/Exp" + str(exp) + "/Raw_data/grid_frame" + str(frame), 'r')
        f.readline()
        f.readline()  # to skip the first 2 header lines

        for line in f:  # each line has information about one of the proteins in the grid
            line = line.strip()
            line_lst = line.split(sep=",")
            x = int(line_lst[0]) * square_pix_size  # resizing x and y locations
            y = int(line_lst[1]) * square_pix_size

            # color selection by complexID
            if line_lst[4] == '0':  # if the protein is not in a complex (3+)
                if pix[x, y] != (255, 255, 255):  # so locations with both a protein and a complex would be colored
                    # with the complex color
                    continue
                color = (228, 221, 222)
            else:  # if the protein is in a complex
                idx = (int(line_lst[4]) - 1) % 10
                color = color_lst[idx]

            # paint the square in the selected color
            for i in range(square_pix_size):
                for j in range(square_pix_size):
                    pix[x + i, y + j] = color

        for i in range(0, size, square_pix_size):  # to create white spaces between each square in the grid
            for j in range(size):
                pix[i, j] = (255, 255, 255)
                pix[j, i] = (255, 255, 255)
        f.close()

        # correct the image format so it will be ready for mp4 format
        img = img.resize((1000, 1000))
        img = img.convert("RGBA")

        # save the image
        img.save(path + "/Exp" + str(exp) + "/Raw_data/frame" + str(frame) + ".png")

# create a visualization of the grid for each of the frames by isTrans (later those images will turn into a
# movie) #
def vis_by_trans(exp, frames, grid_len, path):
    square_pix_size = 10  # one square in the grid will be square_pix_size x square_pix_size pixels
    size = (square_pix_size * grid_len) + 1  # calculate the image size

    for frame in range(frames + 1):  # create an image for each frame based on the grid data in Raw_data folder -
        # only to exp that was marked for visualization.
        img = Image.new('RGB', (size, size), color=(255, 255, 255))
        pix = img.load()
        f = open(path + "/Exp" + str(exp) + "/Raw_data/grid_frame" + str(frame), 'r')
        f.readline()
        f.readline()  # to skip the first 2 header lines

        for line in f:  # each line has information about one of the proteins in the grid
            line = line.strip()
            line_lst = line.split(sep=",")
            x = int(line_lst[0]) * square_pix_size  # resizing x and y locations
            y = int(line_lst[1]) * square_pix_size

            # color selection by isTrans
            if line_lst[3] == 'True':  # if the protein is in trans interaction
                color = (240, 113, 103)
            elif line_lst[3] == 'False':  # if the protein is not in trans interaction
                color = (228, 221, 222)
            else:
                raise IOError("values read from the Raw Data folder are incorrect")

            # paint the square in the selected color
            for i in range(square_pix_size):
                for j in range(square_pix_size):
                    pix[x + i, y + j] = color

        for i in range(0, size, square_pix_size):  # to create white spaces between each square in the grid
            for j in range(size):
                pix[i, j] = (255, 255, 255)
                pix[j, i] = (255, 255, 255)
        f.close()

        # correct the image format so it will be ready for mp4 format
        img = img.resize((1000, 1000))
        img = img.convert("RGBA")

        # save the image
        img.save(path + "/Exp" + str(exp) + "/Raw_data/frame" + str(frame) + ".png")

# create a visualization of the grid for each of the frames by membrane affiliation (later those images will turn into a
# movie) #
def vis_by_membrane(exp, frames, grid_len, path):
    square_pix_size = 10  # one square in the grid will be square_pix_size x square_pix_size pixels
    size = (square_pix_size * grid_len) + 1  # calculate the image size
    #color_lst = [(100, 181, 246), (241, 120, 112), (27, 38, 59)]  # to translate the colour coded words to RGB colours
    color_lst = [(158, 208, 230), (248, 149, 119), (27, 38, 59)]
    for frame in range(frames + 1):  # create an image for each frame based on the grid data in Raw_data folder -
        # only to exp that was marked for visualization.
        img = Image.new('RGB', (size, size), color=(255, 255, 255))
        pix = img.load()
        f = open(path + "/Exp" + str(exp) + "/Raw_data/grid_frame" + str(frame), 'r')
        f.readline()
        f.readline()  # to skip the first 2 header lines
        for line in f:
            line = line.strip()
            line_lst = line.split(sep=",")
            color = color_lst[int(line_lst[6])]
            x = int(line_lst[0]) * square_pix_size
            y = int(line_lst[1]) * square_pix_size
            for i in range(square_pix_size):
                for j in range(square_pix_size):
                    if pix[x + i, y + j] != (255, 255, 255):
                        pix[x + i, y + j] = color_lst[2]
                    else:
                        pix[x + i, y + j] = color
        for i in range(0, size, square_pix_size):  # to create white spaces between each square in the grid
            for j in range(size):
                pix[i, j] = (255, 255, 255)
                pix[j, i] = (255, 255, 255)
        f.close()

        # correct the image format so it will be ready for mp4 format
        img = img.resize((1000, 1000))
        img = img.convert("RGBA")

        # save the image
        img.save(path + "/Exp" + str(exp) + "/Raw_data/frame" + str(frame) + ".png")

# create a movie from images of visualization #
def mve_from_img(exp, frames, path):
    files = [path + "/Exp"+str(exp) + "/Raw_data/frame" + str(i) + ".png" for i in range(frames + 1)]
    mve = ImageSequenceClip(files, fps=5)
    mve.write_videofile(path+"/Exp"+str(exp)+"/exp"+str(exp)+"_movie.mp4", fps=40)


# analyse 'stst_exp' files to 3 basic graphs - number of cis, number of trans and complex length #
def graph_stat(ExpNum, frames, path):
    cis_amount = pd.DataFrame(columns=['frame', 'cis'])  # data-frame for cis amount graph
    trans_amount = pd.DataFrame(columns=['frame', 'trans'])  # data-frame for trans amount graph
    complex_len = pd.DataFrame(columns=['frame', 'len', 'type'])  # data-frame for avg and max length graph

    for exp in range(1, ExpNum + 1):  # export the data from each of the samples in the exp into the data-frames
        filename = path+"/Exp"+str(exp)+"/Stats_Exp"+str(exp)
        file = open(filename, 'r')
        file.readline()
        file.readline()  # to skip the first 2 lines of header
        index = 1  # to remember the frame number (each line we read is a frame number)
        for line in file:
            line = line.strip()
            line_lst = line.split(sep=',')
            cis_amount = cis_amount.append({'frame': index, 'cis': int(float(line_lst[0]))},
                                           ignore_index=True)
            trans_amount = trans_amount.append({'frame': index, 'trans': int(float(line_lst[1]))},
                                               ignore_index=True)
            complex_len = complex_len.append({'frame': index, 'len': float(line_lst[2]), 'type': 'avg'},
                                             ignore_index=True)
            complex_len = complex_len.append({'frame': index, 'len': int(float(line_lst[3])), 'type': 'max'},
                                             ignore_index=True)
            index += 1
        file.close()

    # setting stuff for seaborn
    sns.set_color_codes('pastel')
    hue_len_col = {'max': 'g', 'avg': 'm'}  # for the complex length graph
    sns.set(style='ticks', rc={'lines.linewidth': 0.7})

    # creating the cis graph
    cis_graph = sns.catplot(x='frame', y='cis', color='b', data=cis_amount, kind='point', height=5, aspect=4)
    cis_graph.fig.suptitle("Number of Cis Interactions per Frame", fontsize=20)
    plt.xlabel("Frame Number", fontsize=16)
    plt.ylabel("Number Of Cis Interactions", fontsize=16)
    plt.savefig(path+'/cis_graph')

    # creating the trans graph
    trans_graph = sns.catplot(x='frame', y='trans', color='r', data=trans_amount, kind='point', height=5, aspect=4)
    trans_graph.fig.suptitle("Number of Trans Interactions per Frame", fontsize=20)
    plt.xlabel("Frame Number", fontsize=16)
    plt.ylabel("Number Of Trans Interactions", fontsize=16)
    plt.savefig(path+'/trans_graph')

    # creating the complex length graph
    plt.clf()
    plt.figure(figsize=(25, 5))  # since this graph is long we will set the image length
    comp_graph = sns.boxplot(x='frame', y='len', hue='type', palette=hue_len_col, data=complex_len, showfliers=False,
                             linewidth=0.6, width=0.9)
    comp_graph.set_title("Average and Maximum Length of Complexes per Frame", fontsize=20)
    plt.xlabel("Frame Number", fontsize=16)
    plt.ylabel("Length of Complexes", fontsize=16)
    plt.savefig(path+'/complex_boxplot')


# create 4 heatmaps of different Kas (cis number, trans number, avg complex len and max complex len).
# In addition create 2 files of database - one with all the information of the 81 heatmap cells and one with only the
# average value of the 81 cells #
def Ka_heatmap(conc, trap, expnum, frames, iso_idn, file_path, Annot_HeatMap, jump):
    if not os.path.isfile(str(file_path) + "/AvgDataBase_toR.csv"):  # if the data bases are not available, create them
        # creating a data base with all the data of 81 exps and their samples
        db = pd.DataFrame(columns=['ka_cis', 'ka_trans', 'exp', 'frame', 'cis_num', 'trans_num', 'avg_len', 'max_len',
                                   'conc', 'iso_idn', 'trap_size'])

        # create a file that the data will later append to. the data will append after every 500 lines to save memory
        db.to_csv(str(file_path) + "/FullDataBase_toR.csv", mode='w', index=False, header=True)

        flag_500 = 0
        for cis in range(0, 9, jump):
            for trans in range(0, 9, jump):
                for exp in range(1, expnum + 1):
                    file_name = file_path + "/Exp_" + str(conc) + "-" + str(trap) + "_CIS" + str(cis) + "_TRANS" + \
                                str(trans) + "_Frames" + str(frames) + "_ExpNo" + str(expnum) + "/Exp" + str(exp) + \
                                "/Stats_Exp" + str(exp)
                    f = open(file_name, 'r')
                    f.readline()
                    f.readline()  # to skip the first 2 lines of header

                    frame = 1  # to keep track on the frame number (every line in the file is a frame)
                    for line in f:
                        if flag_500 >= 500:  # export data to the csv file and re-initiate the dataframe ('db' var)
                            db.to_csv(str(file_path) + "/FullDataBase_toR.csv", mode='a', index=False, header=False)
                            flag_500 = 0
                            db = db[0:0]

                        flag_500 += 1
                        line = line.strip()
                        line_lst = line.split(sep=",")
                        db = db.append({'ka_cis': cis, 'ka_trans': trans, 'exp': exp, 'frame': frame,
                                        'cis_num': int(float(line_lst[0])), 'trans_num': int(float(line_lst[1])),
                                        'avg_len': (float(line_lst[2])), 'max_len': int(float(line_lst[3])),
                                        'conc': conc, 'iso_idn': str(iso_idn), 'trap_size': trap},
                                       ignore_index=True)
                        frame += 1
        db.to_csv(str(file_path) + "/FullDataBase_toR.csv", mode='a', index=False, header=False)

        db = pd.read_csv(str(file_path) + "/FullDataBase_toR.csv", usecols=['ka_cis', 'ka_trans', 'exp', 'frame',
                                                                            'cis_num', 'trans_num', 'avg_len',
                                                                            'max_len'])

        # create a summary db - with the average value of the samples in each exp and frame
        sum_db = pd.DataFrame(columns=['ka_cis', 'ka_trans', 'frame', 'avg_cis_num', 'avg_trans_num', 'avg_avg_len',
                                       'avg_max_len', 'conc', 'iso_idn', 'trap_size'])
        for cis in range(0, 9, jump):
            for trans in range(0, 9, jump):
                for frame in range(1, frames + 1):
                    cond = (db['ka_cis'] == cis) & (db['ka_trans'] == trans) & (db['frame'] == frame)
                    avg_cis_num = db[cond]['cis_num'].mean()
                    avg_trans_num = db[cond]['trans_num'].mean()
                    avg_avg_len = db[cond]['avg_len'].mean()
                    avg_max_len = db[cond]['max_len'].mean()
                    sum_db = sum_db.append({'ka_cis': cis, 'ka_trans': trans, 'frame': frame, 'avg_cis_num': avg_cis_num,
                                            'avg_trans_num': avg_trans_num, 'avg_avg_len': avg_avg_len,
                                            'avg_max_len': avg_max_len, 'conc': conc, 'iso_idn': str(iso_idn),
                                            'trap_size': trap}, ignore_index=True)
        sum_db.to_csv(file_path + "/AvgDataBase_toR.csv")  # save the data as csv file for easy access next time
        db = 0  # erase the full 'db' var from memory
        print("Finished creating the data bases")

        # creating directories to save the images
        os.mkdir(file_path + "/Heatmap_Cis")
        os.mkdir(file_path + "/Heatmap_Trans")
        os.mkdir(file_path + "/Heatmap_AvgLen")
        os.mkdir(file_path + "/Heatmap_MaxLen")

    else:
        print("Using existing databases....")
        sum_db = pd.read_csv(file_path + "/AvgDataBase_toR.csv")

    # create the heatmaps for each frame
    for frame in range(1, frames + 1):
        tmp_db = sum_db[sum_db['frame'] == frame]  # create a sub-db of the frame information
        heatmap_typs = ['avg_cis_num', 'avg_trans_num', 'avg_avg_len', 'avg_max_len']
        vmax_dict = {'avg_cis_num': 100, 'avg_trans_num': 100, 'avg_avg_len': 25, 'avg_max_len': 30}  # dictionary
        # with values for the heatmap high-value

        for typ in heatmap_typs:  # create 4 heatmaps for every frame
            pivot_tmp_db = tmp_db.pivot(index='ka_cis', columns='ka_trans', values=typ)
            ax = sns.heatmap(pivot_tmp_db, square=True, linewidths=0.6, cbar=True, annot=Annot_HeatMap, cmap="OrRd",
                             vmin=0, vmax=vmax_dict[typ])
            ax.invert_yaxis()
            if typ == heatmap_typs[0]:
                plt.xlabel(r'$\Delta{G}_D$' + '-Trans', fontsize=14)
                plt.ylabel(r'$\Delta{G}_D$' + '-Cis', fontsize=14)
                plt.title("Average Number of Cis Interactions")
                plt.savefig(file_path + "/Heatmap_Cis/frame" + str(frame) + ".png")
            elif typ == heatmap_typs[1]:
                plt.xlabel(r'$\Delta{G}_D$' + '-Trans', fontsize=14)
                plt.ylabel(r'$\Delta{G}_D$' + '-Cis', fontsize=14)
                plt.title("Average Number of Trans Interactions")
                plt.savefig(file_path + "/Heatmap_Trans/frame" + str(frame) + ".png")
            elif typ == heatmap_typs[2]:
                plt.xlabel(r'$\Delta{G}_D$' + '-Trans', fontsize=14)
                plt.ylabel(r'$\Delta{G}_D$' + '-Cis', fontsize=14)
                plt.title("Average Complex length")
                plt.savefig(file_path + "/Heatmap_AvgLen/frame" + str(frame) + ".png")
            elif typ == heatmap_typs[3]:
                plt.xlabel(r'$\Delta{G}_D$' + '-Trans', fontsize=14)
                plt.ylabel(r'$\Delta{G}_D$' + '-Cis', fontsize=14)
                plt.title("Average Maximum Complex length")
                plt.savefig(file_path + "/Heatmap_MaxLen/frame" + str(frame) + ".png")
            plt.close()


# create a line plots in changing Ka (cis/tran) while the other Ka is constant. this function uses the data-bases
# created by the 'Ka_heatmap' function #
def Ka_lineplot(const_cis, const_trans, file_path):
    db = pd.read_csv(file_path + "/AvgDataBase_toR.csv")  # read the data base created with the 'KD_heatmap' function

    # setting stuff for seaborn
    sns.set_theme(style="ticks")
    palette = sns.color_palette("rocket_r")

    # creating graphs with constant cis and changing trans
    tmp_db = db[db['ka_cis'] == const_cis]
    fig1 = sns.relplot(data=tmp_db, x="frame", y="avg_cis_num", hue="ka_trans", kind='line')
    plt.savefig(file_path + "/NumberOfCisPerFrame_ConstCis" + str(const_cis))
    fig2 = sns.relplot(data=tmp_db, x="frame", y="avg_trans_num", hue="ka_trans", kind='line')
    plt.savefig(file_path + "/NumberOfCisPerFrame_ConstCis" + str(const_cis))

    # creating graphs with constant trans and changing cis
    tmp_db = db[db['ka_trans'] == const_trans]
    fig1 = sns.relplot(data=tmp_db, x="frame", y="avg_cis_num", hue="ka_cis", kind='line')
    plt.savefig(file_path + "/NumberOfCisPerFrame_ConstTrans" + str(const_trans))
    fig2 = sns.relplot(data=tmp_db, x="frame", y="avg_trans_num", hue="ka_cis", kind='line')
    plt.savefig(file_path + "/NumberOfCisPerFrame_ConstTrans" + str(const_trans))


# create a database for R analysis #
def Histogram_db_creator(conc, trap, expnum, frames, file_path):
    db = pd.DataFrame(columns=['Exp', 'Zipper_len', 'Ka'])
    kas = [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 6], [7, 7], [8, 8]]  # which cis and trans combinations to
    # export. need to be a list of lists: [[cis,trans], [cis,trans], ...]
    for pair in kas:
        for exp in range(1, expnum + 1):
            file_name = file_path + "/Exp_" + str(conc) + "-" + str(trap) + "_CIS" + str(pair[0]) + "_TRANS" + str(
                pair[1]) + "_Frames" + str(frames) + "_ExpNo" + str(expnum) + "/Exp" + str(
                exp) + "/Full_Complex_Length_Exp" + str(exp)
            f = open(file_name, 'r')
            for position, line in enumerate(f):
                if position == frames - 1:
                    line = line.strip()
                    line_lst = line.split(sep=",")
                    line_lst.pop(-1)
                    for i in range(len(line_lst)):
                        db = db.append({'Exp': exp, 'Zipper_len': int(line_lst[i]), 'Ka': str(pair[0]) + str(pair[1])}, ignore_index=True)
            f.close()

    db.to_csv(file_path + '/data_for_hist.csv')
