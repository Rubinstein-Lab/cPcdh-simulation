import numpy as np
import random
import math
import pandas as pd


# grid axis directions #
#               y
#      |- - - - ->
#      |
#      |
#    x V

# an empty cell will have a value of 0 #
# protein without cis/tran will have a value of -1 in the respected field #

########################################################################################################################


# class definitions and information #

class Complex:
    def __init__(self, prot_lst, directions):
        # input check #
        if not (isinstance(directions, list)):
            raise ValueError("trying to create complex with wrong directions")
        if len(directions) != 2:
            raise ValueError("trying to create complex with wrong directions")
        if not (isinstance(prot_lst, list)):
            raise ValueError("trying to create complex with wrong proteins")

        # create self #
        self.direction = directions  # a list with 2 direction values, representing the directions of the proteins in
        # the complex
        self.proteins = prot_lst  # list containing all the proteins within the complex
        self.col = 0  # for visualization through the 'save_grid' function

    def __len__(self):  # the length of the complex is defined as the number of proteins within it
        return len(self.proteins)

    def __repr__(self):
        return "complex"


class Protein:
    def __init__(self, location, iso_type, complex, direction, cis=-1):
        # input check #
        if not (isinstance(location, list)):
            raise ValueError("trying to create protein with wrong location")
        if len(location) != 3:
            raise ValueError("trying to create protein with wrong location")
        if location[2] not in [0, 1]:
            raise ValueError("trying to create protein with wrong location")
        if not (isinstance(complex, Complex)):
            raise ValueError("trying to create protein with wrong complex")
        if not (isinstance(direction, str)):
            raise ValueError("trying to create protein with wrong direction")
        if cis != -1:
            if not isinstance(cis, Protein):
                raise ValueError("trying to create protein with wrong cis type")

        # create self #
        self.location = location  # the protein location in the grid as [x,y,membrane]. the membrane is 0/1.
        self.type = iso_type  # the isoform type as an integer
        self.cis = cis  # the protein it interacts with in cis. if no cis - this value will be -1.
        self.trans = -1  # the protein it interacts with in trans. if no trans - this value will be -1.
        self.direction = direction  # the protein direction. one of eight options: 'N'/'E'/'W'/'S'/'SW'/'SE'/'NW'/'NE.
        self.complex = complex  # the complex this protein is apart of.
        self.flag = False  # this field is used during analysis in 'analyse_grid' to make sure each protein will only
        # be checked once.

        if location[2] == 0:  # can be used for the debug process
            self.num = 1
        if location[2] == 1:  # can be used for the debug process
            self.num = 4

    def __repr__(self):  # name of the protein in the grid is the direction - can be modified for debug processes
        return self.direction

    def __add__(self, other):  # can be used for the debug process
        if isinstance(other, int):
            return self.num + other
        if isinstance(other, Protein):
            return self.num + other.num

    def __radd__(self, other):  # can be used for the debug process
        if isinstance(other, int):
            return self.num + other
        if isinstance(other, Protein):
            return self.num + other.num


########################################################################################################################

# general functions #


# corrects the location x or(!) y value in case it's going out of the grid range #
def adjust_xy(val, grid_len):
    if val < 0:
        val = grid_len + val
    if val >= grid_len:
        val = val - grid_len
    return val


# return the inverse direction of a given direction
def inverse_direction(direction):
    dic = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E', 'NW': 'SE', 'NE': 'SW', 'SE': 'NW', 'SW': 'NE'}
    inve_dir = dic[direction]
    return inve_dir


# calculate the diffusion trap locations #
def dif_trap_loc(grid_len,
                 dif_trap):  # create a list with [xy_min, xy_max] that reprasnt the start and end of the trap x,y
    trap_len = math.floor(math.sqrt(dif_trap * (grid_len ** 2)))
    xy_min = (grid_len - trap_len) / 2
    xy_max = grid_len - xy_min
    xy_min = math.floor(xy_min)
    xy_max = math.floor(xy_max)
    if not isinstance(xy_min, int) or not isinstance(xy_max, int):
        raise SystemError("could not find the trap locations")
    trap_locations = [xy_min, xy_max]
    return trap_locations


# perform a non-uniform steps per frame function #
def StepNonUni(frame):
    if frame < 25:
        return 1000
    elif frame < 50:
        return 5000
    elif frame < 75:
        return 100000
    else:
        return 300000


# look for proteins that have trans interactions next to a specific cell, used for 'analyse_grid' function #
def find_around(x, y, grid):
    grid_len = len(grid)
    around_lst = list()
    for up in range(-1, 2):
        for down in range(-1, 2):
            if up == 0 and down == 0:
                continue
            new_x = adjust_xy(x + up, grid_len)
            new_y = adjust_xy(y + down, grid_len)
            # this lines check only one membrane but look for trans interaction, make sure we won't therefore true for
            # both check the same proteins more then once through the 'flag' field.
            if grid[new_x, new_y, 0] != 0 and grid[new_x, new_y, 0].trans != -1 and grid[new_x, new_y, 0].flag == False:
                around_lst.append([new_x, new_y])
    return around_lst


# analyse and find clusters of trans interactions for the 'analyse_grid' function #
def cluster(x, y, grid):
    cluster_lst = [[x, y]]
    around_lst = find_around(x, y, grid)
    grid[x, y, 0].flag = True
    if len(around_lst) > 1:  # if more than 1 protein is located near our query
        for ar in around_lst:
            grid[ar[0], ar[1], 0].flag = True  # to flag checked cells
        cluster_lst = cluster_lst + around_lst
        for near in around_lst:
            next_around = cluster(near[0], near[1], grid)
            cluster_lst = cluster_lst + next_around
    return cluster_lst


########################################################################################################################

# functions for the initiation process #


# create an empty 3D array in the size of (len x len x 2) #
def create_empty_grid(grid_len):  # create a grid with zeros
    grid = np.zeros((grid_len, grid_len, 2), dtype=Protein)
    return grid


# locates empty spots for cis-dimers insertion during initial grid formation #
def find_spot_for_dimer(grid, membrane):
    grid_len = len(grid)
    x1 = random.choice(range(grid_len))  # choose random x&y values
    y1 = random.choice(range(grid_len))
    while grid[x1, y1, membrane] != 0:  # keep choosing randomly until an empty spot is located for one(!) of the two
        # proteins in the cis-dimer
        x1 = random.choice(range(grid_len))
        y1 = random.choice(range(grid_len))
    direction = random.choice([['N', 'S'], ['W', 'E'], ['SE', 'NW'], ['SW', 'NE']])  # choose a random direction for
    # the cis-dimer

    # based of the first location and the cis-dimer find if the corresponding location of the second protein is empty
    if direction == ['N', 'S']:
        x2 = adjust_xy(x1 + 1, grid_len)
        y2 = y1
        if grid[x2, y2, membrane] == 0:
            return [[x1, y1, membrane], [x2, y2, membrane]], direction  # if the second location is free the function
            # will return both locations and the cis-dimer direction
        x2 = adjust_xy(x1 - 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)  # using recursion if the second location is not empty to
            # restart the process of random location

    if direction == ['W', 'E']:
        x2 = x1
        y2 = adjust_xy(y1 + 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        y2 = adjust_xy(y1 - 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)

    if direction == ['SE', 'NW']:
        x2 = adjust_xy(x1 - 1, grid_len)
        y2 = adjust_xy(y1 - 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x2, y1, membrane] != 0:  # to make sure cis-dimers won't intersect
                if grid[x2, y1, membrane].direction in ['NE', 'SW']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        x2 = adjust_xy(x1 + 1, grid_len)
        y2 = adjust_xy(y1 + 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x2, y1, membrane] != 0:  # to make sure cis-dimers won't intersect
                if grid[x2, y1, membrane].direction in ['NE', 'SW']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)

    if direction == ['SW', 'NE']:
        x2 = adjust_xy(x1 - 1, grid_len)
        y2 = adjust_xy(y1 + 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x2, y1, membrane] != 0:  # to make sure cis-dimers won't intersect
                if grid[x2, y1, membrane].direction in ['NW', 'SE']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        x2 = adjust_xy(x1 + 1, grid_len)
        y2 = adjust_xy(y1 - 1, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x2, y1, membrane] != 0:  # to make sure cis-dimers won't intersect
                if grid[x2, y1, membrane].direction in ['NW', 'SE']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)

    # raise an error if an insertion has been skipped
    raise SystemError("dimer skipped from insertion to grid")


# choose protein isoforms to the dimer #
def choose_types(iso_lst, iso_dict, Cis_Mrx):
    typ1 = random.choice(iso_lst)  # randomly assigning the proteins an isoform type
    while iso_dict[typ1] == 0:  # checking the isoform dictionary inventory to make sure all isoforms
        # will have the input concentration
        typ1 = random.choice(iso_lst)
    typ2 = random.choice(iso_lst)
    while iso_dict[typ2] == 0:
        typ2 = random.choice(iso_lst)

    # try creating a cis dimer according to the affinity
    # identify the cis affinity from the affinity matrix
    if not pd.isnull(Cis_Mrx.iloc[(typ1 - 1), (typ2 - 1)]):
        KA_cis = Cis_Mrx.iloc[(typ1 - 1), (typ2 - 1)]
    else:
        KA_cis = Cis_Mrx.iloc[(typ2 - 1), (typ1 - 1)]

    deltaG = np.exp(KA_cis) / (np.exp(KA_cis) + 1)  # apply the affinity function
    rnd = random.random()
    if rnd > deltaG:  # if we decide not to create this interaction then restart
        return None, None  # else we will restart the grid initiation process

    return typ1, typ2


# insert a protein to a grid location and return the new grid #
def insert_protein_to_grid(grid, prot):
    grid[prot.location[0], prot.location[1], prot.location[2]] = prot
    return grid


# the main initiation function - this is the only function that the 'Initiator.py' excess #
# creates a randomly chosen grid of cis-dimers based on the used inputs #
def init_grid(grid_len, isoAmount, num_SameIso, num_DifIso, num_TotIso, Cis_Mrx):
    # creates an empty grid
    grid = create_empty_grid(grid_len)

    # creates an isoform list (represented as integers) to each of the membranes based on the used input
    membrane0_iso_lst = list(range(1, num_SameIso + 1)) + list(range(num_SameIso + 1, num_SameIso + num_DifIso + 1))
    membrane1_iso_lst = list(range(1, num_SameIso + 1)) + list(range(num_SameIso + num_DifIso + 1,
                                                                     num_SameIso + (num_DifIso * 2) + 1))

    # transform the lists into a dictionaries with amounts. this dictionary will be updated after every insertion and
    # allow to make sure all the isoforms has been inserted at the correct concentration
    membrane0_iso_dict = dict()
    for i in membrane0_iso_lst:
        membrane0_iso_dict[i] = isoAmount
    membrane1_iso_dict = dict()
    for i in membrane1_iso_lst:
        membrane1_iso_dict[i] = isoAmount

    # create unique ID numbers for the proteins and the complexes
    all_protein_id = range(1, (num_TotIso * isoAmount * 2) + 1)
    all_protein = ['prot' + str(i) for i in all_protein_id]  # list of all the proteins names
    all_complex_id = range(1, (num_TotIso * isoAmount) + 1)
    all_complex = ['comp' + str(i) for i in all_complex_id]  # list of all the cis-dimers names

    prot_index = 0  # for keeping the protein ID name index as creating the proteins

    # start creating cis-dimers
    # inserting cis-dimers to membrane 0
    for comp_index in range((len(all_complex)) // 2):  # first half of complexes ID are assigned to membrane 0
        comp_name = all_complex[comp_index]
        loc = find_spot_for_dimer(grid, 0)  # finding a location for the cis-dimer
        prot1_name = all_protein[prot_index]
        prot2_name = all_protein[prot_index + 1]

        # assigning isoforms types to the proteins in cis-dimer
        typ1 = None
        flag = 0
        while typ1 is None:
            typ1, typ2 = choose_types(membrane0_iso_lst, membrane0_iso_dict, Cis_Mrx)
            flag += 1  # to check if we are stuck on the loop
            if flag >= 10000:
                prot_left = 0  # count how many proteins hasn't found a cis-dimer yet
                flag_low = True
                for inx in membrane0_iso_lst:
                    if membrane0_iso_dict[inx] >= (isoAmount * 0.1):
                        print(inx)
                        flag_low = False
                    prot_left = prot_left + membrane0_iso_dict[inx]
                if flag_low:
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print(prot_left)

        print(typ1, typ2)

        membrane0_iso_dict[typ1] -= 1  # updating the isoform dictionary inventory
        membrane0_iso_dict[typ2] -= 1

        # creating the complex and proteins - with the classes defined before (Complex & Protein)
        locals()[comp_name] = Complex([0, 0], loc[1])
        locals()[prot1_name] = Protein(loc[0][0], typ1, locals()[comp_name], loc[1][0])
        locals()[prot2_name] = Protein(loc[0][1], typ2, locals()[comp_name], loc[1][1])
        locals()[comp_name].proteins = [locals()[prot1_name], locals()[prot2_name]]  # updating the complex proteins
        # list field of the class
        locals()[prot1_name].cis = locals()[prot2_name]  # updating the protein's cis interactions field of the class
        locals()[prot2_name].cis = locals()[prot1_name]

        # adding the proteins to the grid
        grid = insert_protein_to_grid(grid, locals()[prot1_name])
        grid = insert_protein_to_grid(grid, locals()[prot2_name])

        # updating the 'all_protein' and 'all_complex' lists - those lists will use later for accessing the proteins
        # and the complexes
        all_protein[prot_index] = locals()[prot1_name]
        all_protein[prot_index + 1] = locals()[prot2_name]
        all_complex[comp_index] = locals()[comp_name]

        prot_index += 2  # updating the protein ID name index after every cis-dimer insertion

    # inserting cis-dimers to membrane 1 - similar
    for comp_index in range((len(all_complex)) // 2, len(all_complex)):  # second half of complexes ID are assigned to
        # membrane 1
        comp_name = all_complex[comp_index]
        loc = find_spot_for_dimer(grid, 1)
        prot1_name = all_protein[prot_index]
        prot2_name = all_protein[prot_index + 1]

        # assigning isoforms types to the proteins in cis-dimer
        typ1 = None
        flag = 0
        while typ1 is None:
            typ1, typ2 = choose_types(membrane1_iso_lst, membrane1_iso_dict, Cis_Mrx)
            flag += 1  # to check if we are stuck on the loop
            if flag >= 10000:
                prot_left = 0  # count how many proteins hasn't found a cis-dimer yet
                flag_low = True
                for inx in membrane1_iso_lst:
                    if membrane1_iso_dict[inx] >= (isoAmount * 0.1):
                        print(inx)
                        flag_low = False
                    prot_left = prot_left + membrane1_iso_dict[inx]
                if flag_low:
                    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                print(prot_left)

        print(typ1, typ2)

        membrane1_iso_dict[typ1] -= 1  # updating the isoform dictionary inventory
        membrane1_iso_dict[typ2] -= 1

        locals()[comp_name] = Complex([0, 0], loc[1])
        locals()[prot1_name] = Protein(loc[0][0], typ1, locals()[comp_name], loc[1][0])
        locals()[prot2_name] = Protein(loc[0][1], typ2, locals()[comp_name], loc[1][1])
        locals()[comp_name].proteins = [locals()[prot1_name], locals()[prot2_name]]
        locals()[prot1_name].cis = locals()[prot2_name]
        locals()[prot2_name].cis = locals()[prot1_name]

        grid = insert_protein_to_grid(grid, locals()[prot1_name])
        grid = insert_protein_to_grid(grid, locals()[prot2_name])

        all_protein[prot_index] = locals()[prot1_name]
        all_protein[prot_index + 1] = locals()[prot2_name]
        all_complex[comp_index] = locals()[comp_name]

        prot_index += 2

    # after the initiation process we will get an initial random grid (the array with the proteins inside),
    # and a 2 lists that include all the proteins and all the complexes
    return grid, all_protein, all_complex


# locates empty spots for cis-dimers insertion during initial grid formation #
def find_spot_for_monomer(grid, membrane):
    grid_len = len(grid)
    x1 = random.choice(range(grid_len))  # choose random x&y values
    y1 = random.choice(range(grid_len))
    while grid[x1, y1, membrane] != 0:  # keep choosing randomly until an empty spot is located for the monomer
        x1 = random.choice(range(grid_len))
        y1 = random.choice(range(grid_len))
    direction = random.choice(['N', 'S', 'E', 'W', 'NE', 'NW', 'SE', 'SW'])  # choose a random direction for
    # the monomer

    return [x1, y1, membrane], direction


# the main initiation function - this is the only function that the 'Initiator.py' excess #
# creates a randomly chosen grid of cis-dimers and monomers based on the used inputs #
def init_grid2(grid_len, isoAmount, num_SameIso, num_DifIso, num_TotIso, Cis_Mrx, trap_locations):
    # creates an empty grid
    grid = create_empty_grid(grid_len)

    # creates an isoform list (represented as integers) to each of the membranes based on the used input
    membrane0_iso_lst = list(range(1, num_SameIso + 1)) + list(range(num_SameIso + 1, num_SameIso + num_DifIso + 1))
    membrane1_iso_lst = list(range(1, num_SameIso + 1)) + list(range(num_SameIso + num_DifIso + 1,
                                                                     num_SameIso + (num_DifIso * 2) + 1))

    # transform the lists into a dictionaries with amounts. this dictionary will be updated after every insertion and
    # allow to make sure all the isoforms has been inserted at the correct concentration
    membrane0_iso_dict = dict()
    for i in membrane0_iso_lst:
        membrane0_iso_dict[i] = isoAmount
    membrane1_iso_dict = dict()
    for i in membrane1_iso_lst:
        membrane1_iso_dict[i] = isoAmount

    # create unique ID numbers for the proteins and the complexes
    all_protein_id = range(1, (num_TotIso * isoAmount * 2) + 1)
    all_protein = ['prot' + str(i) for i in all_protein_id]  # list of all the proteins names
    all_complex_id = range(1, (num_TotIso * isoAmount * 2) + 1)
    all_complex = ['comp' + str(i) for i in all_complex_id]  # list of all the cis-dimers names

    prot_index = 0  # for keeping the protein ID name index as creating the proteins

    # start creating monomers
    # inserting monomers to membrane 0
    for comp_index in range((len(all_complex)) // 2):  # first half of complexes ID are assigned to membrane 0
        comp_name = all_complex[comp_index]
        loc1, direction1 = find_spot_for_monomer(grid, 0)  # finding a location for the monomer
        prot1_name = all_protein[prot_index]

        typ1 = random.choice(membrane0_iso_lst)  # randomly assigning the proteins an isoform type
        while membrane0_iso_dict[typ1] == 0:  # checking the isoform dictionary inventory to make sure all isoforms
            # will have the input concentration
            typ1 = random.choice(membrane0_iso_lst)

        membrane0_iso_dict[typ1] -= 1  # updating the isoform dictionary inventory

        # creating the complex and monomer - with the classes defined before (Complex & Protein)
        locals()[comp_name] = Complex([0], [direction1, inverse_direction(direction1)])
        locals()[prot1_name] = Protein(loc1, typ1, locals()[comp_name], direction1)
        locals()[comp_name].proteins = [locals()[prot1_name]]  # updating the complex proteins list field of the class

        # adding the proteins to the grid
        grid = insert_protein_to_grid(grid, locals()[prot1_name])

        # updating the 'all_protein' and 'all_complex' lists - those lists will use later for accessing the proteins
        # and the complexes
        all_protein[prot_index] = locals()[prot1_name]
        all_complex[comp_index] = locals()[comp_name]

        prot_index += 1  # updating the protein ID name index after every monomer insertion

    # inserting cis-dimers to membrane 1 - similar
    for comp_index in range((len(all_complex)) // 2, len(all_complex)):  # second half of complexes ID are assigned to
        # membrane 1
        comp_name = all_complex[comp_index]
        loc1, direction1 = find_spot_for_monomer(grid, 1)  # finding a location for the monomer
        prot1_name = all_protein[prot_index]

        typ1 = random.choice(membrane1_iso_lst)  # randomly assigning the proteins an isoform type
        while membrane1_iso_dict[typ1] == 0:  # checking the isoform dictionary inventory to make sure all isoforms
            # will have the input concentration
            typ1 = random.choice(membrane1_iso_lst)

        membrane1_iso_dict[typ1] -= 1  # updating the isoform dictionary inventory

        # creating the complex and monomer - with the classes defined before (Complex & Protein)
        locals()[comp_name] = Complex([0], [direction1, inverse_direction(direction1)])
        locals()[prot1_name] = Protein(loc1, typ1, locals()[comp_name], direction1)
        locals()[comp_name].proteins = [locals()[prot1_name]]  # updating the complex proteins list field of the class

        # adding the proteins to the grid
        grid = insert_protein_to_grid(grid, locals()[prot1_name])

        # updating the 'all_protein' and 'all_complex' lists - those lists will use later for accessing the proteins
        # and the complexes
        all_protein[prot_index] = locals()[prot1_name]
        all_complex[comp_index] = locals()[comp_name]

        prot_index += 1  # updating the protein ID name index after every monomer insertion

    # after all the proteins has been inserted as monomers, the simulation will use random actions in order to form
    # cis-dimers in each of the membranes separately
    for t in range(10000000):
        what_move = random.choice(['move', 'turn', 'create_cis', 'break_cis'])

        if what_move == 'move':
            grid = move(grid, all_protein, trap_locations)
        elif what_move == 'turn':
            grid = turn(grid, all_protein, trap_locations)
        elif what_move == 'create_cis':
            grid, all_complex, all_protein = create_cis(grid, all_protein, all_complex, Cis_Mrx)
        elif what_move == 'break_cis':
            grid, all_complex, all_protein = break_cis(grid, all_protein, all_complex, Cis_Mrx, 00, t)

    # after the initiation process we will get an initial random grid (the array with the proteins inside),
    # and a 2 lists that include all the proteins and all the complexes
    return grid, all_protein, all_complex


########################################################################################################################

# functions for the 6 actions - move, turn, create cis, break cis, create trans, break trans #


# move #

# sorting a proteins list according to the move direction. from first to move to last to move #
# this will allow to move all the proteins in a complex without running each other in the process #
def sort_prot_by_direction(prot_lst, direction):
    index = 0
    if direction in ['N', 'NW', 'NE']:
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[0])
        for i in range(len(sorted_lst) - 1):  # change the start index for complexes that are crossing the edges
            if sorted_lst[i + 1].location[0] - sorted_lst[i].location[0] not in [0, 1]:
                index = i + 1
    elif direction in ['S', 'SW', 'SE']:
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[0], reverse=True)
        for i in range(len(sorted_lst) - 1):  # change the start index for complexes that are crossing the edges
            if sorted_lst[i].location[0] - sorted_lst[i + 1].location[0] not in [0, 1]:
                index = i + 1
    elif direction == 'E':
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[1], reverse=True)
        for i in range(len(sorted_lst) - 1):  # change the start index for complexes that are crossing the edges
            if sorted_lst[i].location[1] - sorted_lst[i + 1].location[1] not in [0, 1]:
                index = i + 1
    elif direction == 'W':
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[1])
        for i in range(len(sorted_lst) - 1):  # change the start index for complexes that are crossing the edges
            if sorted_lst[i + 1].location[1] - sorted_lst[i].location[1] not in [0, 1]:
                index = i + 1

    # this function return the sorted proteins list and the starting index in the list
    return sorted_lst, index


# check if one protein can move in a given direction #
# same_dir should be 'True' if we move to a direction that is in the direction of the complex #
def check_one_move(grid, prot, direction, same_dir=False):
    grid_len = len(grid)
    side_to_loc = {'N': [-1, 0], 'S': [1, 0], 'W': [0, -1], 'E': [0, 1], 'NE': [-1, 1], 'NW': [-1, -1], 'SE': [1, 1],
                   'SW': [1, -1]}  # gives the coordination different from a cell to a nearby cell in that direction
    dif_x = side_to_loc[direction][0]
    dif_y = side_to_loc[direction][1]
    x_check = adjust_xy(prot.location[0] + dif_x, grid_len)  # calculate the x,y of the cell we need to check
    y_check = adjust_xy(prot.location[1] + dif_y, grid_len)

    if grid[x_check, y_check, prot.location[2]] == 0:  # if the cell we are checking is empty, we need to further
        # check for a possible intersection between two complexes

        if not same_dir:
            if prot.direction in ['N', 'E', 'S', 'W']:
                return True
            else:
                if prot.location[2] == 0:
                    if grid[prot.location[0], prot.location[1], 1] == 0:  # check if the other membrane has a complex
                        return True
                    elif len(grid[prot.location[0], prot.location[1], 1].complex) < 3:
                        return True
                else:
                    if grid[prot.location[0], prot.location[1], 0] == 0:
                        return True
                    elif len(grid[prot.location[0], prot.location[1], 0].complex) < 3:
                        return True

        if same_dir:
            if direction in ['N', 'E', 'S', 'W']:
                return True
            else:
                if grid[x_check, prot.location[1], prot.location[2]] == 0:  # check that both sides are empty as well
                    if grid[prot.location[0], y_check, prot.location[2]] == 0:
                        return True

    # if one of the above statements is false - the move cannot happen
    return False


# move one protein in the grid in a specific direction #
def move_one(grid, prot, direction):
    grid_len = len(grid)
    side_to_loc = {'N': [-1, 0], 'S': [1, 0], 'W': [0, -1], 'E': [0, 1], 'NE': [-1, 1], 'NW': [-1, -1], 'SE': [1, 1],
                   'SW': [1, -1]}  # gives the coordination different from a cell to a nearby cell in that direction
    # calculate the new location
    x_new = adjust_xy(prot.location[0] + side_to_loc[direction][0], grid_len)
    y_new = adjust_xy(prot.location[1] + side_to_loc[direction][1], grid_len)
    z_new = prot.location[2]
    old_loc = prot.location

    # change the grid based on the move and the protein location field
    prot.location = [x_new, y_new, z_new]
    grid[x_new, y_new, z_new] = prot
    grid[old_loc[0], old_loc[1], old_loc[2]] = 0
    return grid


# the main move function - this is the only function that the 'Initiator.py' excess # this function choose a random
# protein and moving direction and test if the move is possible - if possible it execute the move, if not possible
# return the same inputted grid #
def move(grid, all_protein, trap_locations):
    ran_prot = random.choice(all_protein)  # choose a random protein from the proteins list
    complex = ran_prot.complex  # accessing the protein's complex
    if not isinstance(complex, Complex):
        raise SystemError("trying to move a non complex type")
    if len(complex) >= (2 * (len(grid) - 1)):  # if the chosen zipper is in the size of the grid length - he won't move
        return grid
    direction = random.choice(['N', 'S', 'W', 'E', 'NW', 'NE', 'SW', 'SE'])  # choose a random moving direction

    side_to_loc = {'N': [-1, 0], 'S': [1, 0], 'W': [0, -1], 'E': [0, 1], 'NE': [-1, 1], 'NW': [-1, -1], 'SE': [1, 1],
                   'SW': [1, -1]}  # gives the coordination different from a cell to a nearby cell in that direction
    # calculate the new location

    # check if the complex have trans interactions, if so - checks they wont leave the diffusion trap area with this
    # move. if they would leave the diffusion trap with the move it wont execute the move
    for p in complex.proteins:
        if p.trans != -1:
            x_new = adjust_xy(p.location[0] + side_to_loc[direction][0], len(grid))
            y_new = adjust_xy(p.location[1] + side_to_loc[direction][1], len(grid))
            if x_new < trap_locations[0] or x_new > trap_locations[1]:
                return grid
            if y_new < trap_locations[0] or y_new > trap_locations[1]:
                return grid

    # the function separates between two moving options - moving in the direction of the complex and not in the
    # complex direction

    # option one - move in the complex direction
    if direction in complex.direction:
        # for moving in the complex direction there is a significant for checking and moving the proteins in a
        # sequential order. therefore we first sort the protein list of the complex based on the direction
        sorted_prot_lst, index = sort_prot_by_direction(complex.proteins, direction)

        # after sorting, check if the first protein on each membrane can move - if he can, all the other proteins can
        # move as well
        if check_one_move(grid, sorted_prot_lst[index], direction, True):  # check the first protein

            if sorted_prot_lst[index].trans != -1:  # if in trans, check the trans of the first protein
                if check_one_move(grid, sorted_prot_lst[index].trans, direction, True):
                    # if the first trans-dimer can move, start execute the move of the entire complex by order
                    for i in range(index, len(sorted_prot_lst)):
                        grid = move_one(grid, sorted_prot_lst[i], direction)
                    for i in range(0, index):
                        grid = move_one(grid, sorted_prot_lst[i], direction)

            if sorted_prot_lst[index].trans == -1:  # if the first protein not in trans check the trans of the second
                # protein in the sequence
                if sorted_prot_lst[index].cis != -1:
                    if sorted_prot_lst[index].cis.trans != -1:
                        if check_one_move(grid, sorted_prot_lst[index].cis.trans, direction, True):
                            # if the first protein in each membrane can move, start execute the move
                            for i in range(index, len(sorted_prot_lst)):
                                grid = move_one(grid, sorted_prot_lst[i], direction)
                            for i in range(0, index):
                                grid = move_one(grid, sorted_prot_lst[i], direction)

    # option two - move not in the complex direction
    else:
        for prot in complex.proteins:
            if not check_one_move(grid, prot, direction, False):  # if even one protein can't move - don't execute
                return grid
        # only if all proteins can move, start execute the move for all proteins in the complex
        for prot in complex.proteins:
            grid = move_one(grid, prot, direction)

    # return the new grid after the moving process
    return grid


# turn #

# find the expected location of the protein idx(number) cells above a specific protein in a complex # note - this
# function takes into account the complex direction and look for the above in the complex structure and not the grid #
def above_in_chain(grid, prot, idx, mem):
    if prot.direction in ['N', 'S']:
        if mem == 0:
            return grid[adjust_xy(prot.location[0] - idx, len(grid)), prot.location[1], 0]
        else:
            return grid[adjust_xy(prot.location[0] - idx, len(grid)), prot.location[1], 1]

    elif prot.direction in ['NE', 'SW']:
        if mem == 0:
            return grid[adjust_xy(prot.location[0] - idx, len(grid)), adjust_xy(prot.location[1] + idx, len(grid)), 0]
        else:
            return grid[adjust_xy(prot.location[0] - idx, len(grid)), adjust_xy(prot.location[1] + idx, len(grid)), 1]

    elif prot.direction in ['E', 'W']:
        if mem == 0:
            return grid[prot.location[0], adjust_xy(prot.location[1] + idx, len(grid)), 0]
        else:
            return grid[prot.location[0], adjust_xy(prot.location[1] + idx, len(grid)), 1]

    elif prot.direction in ['SE', 'NW']:
        if mem == 0:
            return grid[adjust_xy(prot.location[0] - idx, len(grid)), adjust_xy(prot.location[1] - idx, len(grid)), 0]
        else:
            return grid[adjust_xy(prot.location[0] - idx, len(grid)), adjust_xy(prot.location[1] - idx, len(grid)), 1]


# find the expected location of the protein idx(number) cells below a specific protein in a complex # note - this
# function takes into account the complex direction and look for the below in the complex structure and not the grid #
def below_in_chain(grid, prot, idx, mem):
    if prot.direction in ['N', 'S']:
        if mem == 0:
            return grid[adjust_xy(prot.location[0] + idx, len(grid)), prot.location[1], 0]
        else:
            return grid[adjust_xy(prot.location[0] + idx, len(grid)), prot.location[1], 1]

    elif prot.direction in ['NE', 'SW']:
        if mem == 0:
            return grid[adjust_xy(prot.location[0] + idx, len(grid)), adjust_xy(prot.location[1] - idx, len(grid)), 0]
        else:
            return grid[adjust_xy(prot.location[0] + idx, len(grid)), adjust_xy(prot.location[1] - idx, len(grid)), 1]

    elif prot.direction in ['E', 'W']:
        if mem == 0:
            return grid[prot.location[0], adjust_xy(prot.location[1] - idx, len(grid)), 0]
        else:
            return grid[prot.location[0], adjust_xy(prot.location[1] - idx, len(grid)), 1]

    elif prot.direction in ['SE', 'NW']:
        if mem == 0:
            return grid[adjust_xy(prot.location[0] + idx, len(grid)), adjust_xy(prot.location[1] + idx, len(grid)), 0]
        else:
            return grid[adjust_xy(prot.location[0] + idx, len(grid)), adjust_xy(prot.location[1] + idx, len(grid)), 1]


# check if the complex can move around the turning point in a chosen direction (left or right - randomly chosen) #
def can_complex_turn(grid, prot, trap_locations):
    turn_options = {'N': ['NE', 'NW'], 'S': ['SW', 'SE'], 'E': ['SE', 'NE'], 'W': ['NW', 'SW'], 'SE': ['S', 'E'],
                    'NE': ['N', 'E'], 'SW': ['S', 'W'], 'NW': ['N', 'W']}  # give the new direction based on the old
    # direction and the left/right turn chosen
    check_dif = {'N': [[[0, 1], [0, -1]], [[0, -1], [0, 1]]], 'S': [[[0, 1], [0, -1]], [[0, -1], [0, 1]]],
                 'NE': [[[0, -1], [0, 1]], [[1, 0], [-1, 0]]], 'SW': [[[0, -1], [0, 1]], [[1, 0], [-1, 0]]],
                 'E': [[[1, 0], [-1, 0]], [[-1, 0], [1, 0]]], 'W': [[[1, 0], [-1, 0]], [[-1, 0], [1, 0]]],
                 'SE': [[[0, 1], [0, -1]], [[1, 0], [-1, 0]]], 'NW': [[[0, 1], [0, -1]], [[1, 0], [-1, 0]]]}  # this
    # dictionary will give checking coordination direction for the move based on rnd value, then above/below and then
    # x/y value

    rnd = random.choice([0, 1])  # choose randomly left or right
    grid_len = len(grid)

    # start checking the proteins above the chosen protein
    index = 1  # maintain a count of how many proteins we calculated above the turning point
    while True:
        above_prot_mem0 = above_in_chain(grid, prot, index, 0)
        above_prot_mem1 = above_in_chain(grid, prot, index, 1)
        if above_prot_mem0 == 0 or above_prot_mem0 not in prot.complex.proteins:
            if above_prot_mem1 == 0 or above_prot_mem1 not in prot.complex.proteins:
                break  # if we reached the max above protein - then we tested all above proteins and we can return
                # true for above

        if index > grid_len:  # if too many iterations occurred in the while loop raise an error
            raise SystemError("just a check for now, erase later")

        if above_prot_mem0 in prot.complex.proteins:
            x_check = -1
            y_check = -1
            for i in range(1, index + 1):
                x_check = adjust_xy(above_prot_mem0.location[0] + (i * check_dif[prot.direction][rnd][0][0]),
                                    grid_len)
                y_check = adjust_xy(above_prot_mem0.location[1] + (i * check_dif[prot.direction][rnd][0][1]),
                                    grid_len)
                if grid[x_check, y_check, 0] != 0:
                    return None  # if even one cell is not empty the turn cannot be executed

            # to prevent intersection between complexes we will also check if the other membrane has a complex in the
            # near cell to our new location
            if index > 1 or above_prot_mem0.trans == -1:
                if grid[adjust_xy(x_check - check_dif[prot.direction][rnd][0][0], grid_len),
                        adjust_xy(y_check - check_dif[prot.direction][rnd][0][1], grid_len), 1] != 0:
                    if len((grid[adjust_xy(x_check - check_dif[prot.direction][rnd][0][0], grid_len),
                                 adjust_xy(y_check - check_dif[prot.direction][rnd][0][1], grid_len), 1]).complex) > 2:
                        return None

            # check to see if a trans-dimer leave the diffusion trap area - if so turn cannot be executed
            if above_prot_mem0.trans != -1:
                if x_check < trap_locations[0] or x_check > trap_locations[1]:
                    return None
                elif y_check < trap_locations[0] or y_check > trap_locations[1]:
                    return None

        if above_prot_mem1 in prot.complex.proteins:
            x_check = -1
            y_check = -1
            for i in range(1, index + 1):
                x_check = adjust_xy(above_prot_mem1.location[0] + (i * check_dif[prot.direction][rnd][0][0]),
                                    grid_len)
                y_check = adjust_xy(above_prot_mem1.location[1] + (i * check_dif[prot.direction][rnd][0][1]),
                                    grid_len)
                if grid[x_check, y_check, 1] != 0:
                    return None  # if even one cell is not empty the turn cannot be executed

            # to prevent intersection between complexes we will also check if the other membrane has a complex in the
            # near cell to our new location
            if index > 1 or above_prot_mem1.trans == -1:
                if grid[adjust_xy(x_check - check_dif[prot.direction][rnd][0][0], grid_len),
                        adjust_xy(y_check - check_dif[prot.direction][rnd][0][1], grid_len), 0] != 0:
                    if len((grid[adjust_xy(x_check - check_dif[prot.direction][rnd][0][0], grid_len),
                                 adjust_xy(y_check - check_dif[prot.direction][rnd][0][1], grid_len), 0]).complex) > 2:
                        return None

            # check to see if a trans-dimer leave the diffusion trap area - if so turn cannot be executed
            if above_prot_mem1.trans != -1:
                if x_check < trap_locations[0] or x_check > trap_locations[1]:
                    return None
                elif y_check < trap_locations[0] or y_check > trap_locations[1]:
                    return None
        index += 1

    # start checking the proteins below the chosen protein
    index = 1  # maintain a count of how many proteins we calculated above the turning point
    while True:
        below_prot_mem0 = below_in_chain(grid, prot, index, 0)
        below_prot_mem1 = below_in_chain(grid, prot, index, 1)
        if below_prot_mem0 == 0 or below_prot_mem0 not in prot.complex.proteins:
            if below_prot_mem1 == 0 or below_prot_mem1 not in prot.complex.proteins:
                break  # if we reached the max above protein - then we tested all above proteins and we can return
                # true for above

        if index > grid_len:  # if too many iterations occurred in the while loop raise an error
            raise SystemError("just a check for now, erase later")

        if below_prot_mem0 in prot.complex.proteins:
            x_check = -1
            y_check = -1
            for i in range(1, index + 1):
                x_check = adjust_xy(below_prot_mem0.location[0] + (i * check_dif[prot.direction][rnd][1][0]),
                                    grid_len)
                y_check = adjust_xy(below_prot_mem0.location[1] + (i * check_dif[prot.direction][rnd][1][1]),
                                    grid_len)
                if grid[x_check, y_check, 0] != 0:
                    return None  # if even one cell is not empty the turn cannot be executed

            # to prevent intersection between complexes we will also check if the other membrane has a complex in the
            # near cell to our new location
            if index > 1 or below_prot_mem0.trans == -1:
                if grid[adjust_xy(x_check - check_dif[prot.direction][rnd][1][0], grid_len),
                        adjust_xy(y_check - check_dif[prot.direction][rnd][1][1], grid_len), 1] != 0:
                    if len((grid[adjust_xy(x_check - check_dif[prot.direction][rnd][1][0], grid_len),
                                 adjust_xy(y_check - check_dif[prot.direction][rnd][1][1], grid_len), 1]).complex) > 2:
                        return None

            # check to see if a trans-dimer leave the diffusion trap area - if so turn cannot be executed
            if below_prot_mem0.trans != -1:
                if x_check < trap_locations[0] or x_check > trap_locations[1]:
                    return None
                elif y_check < trap_locations[0] or y_check > trap_locations[1]:
                    return None

        if below_prot_mem1 in prot.complex.proteins:
            x_check = -1
            y_check = -1
            for i in range(1, index + 1):
                x_check = adjust_xy(below_prot_mem1.location[0] + (i * check_dif[prot.direction][rnd][1][0]),
                                    grid_len)
                y_check = adjust_xy(below_prot_mem1.location[1] + (i * check_dif[prot.direction][rnd][1][1]),
                                    grid_len)
                if grid[x_check, y_check, 1] != 0:
                    return None  # if even one cell is not empty the turn cannot be executed

            # to prevent intersection between complexes we will also check if the other membrane has a complex in the
            # near cell to our new location
            if index > 1 or below_prot_mem1.trans == -1:
                if grid[adjust_xy(x_check - check_dif[prot.direction][rnd][1][0], grid_len),
                        adjust_xy(y_check - check_dif[prot.direction][rnd][1][1], grid_len), 0] != 0:
                    if len((grid[adjust_xy(x_check - check_dif[prot.direction][rnd][1][0], grid_len),
                                 adjust_xy(y_check - check_dif[prot.direction][rnd][1][1], grid_len), 0]).complex) > 2:
                        return None

            # check to see if a trans-dimer leave the diffusion trap area - if so turn cannot be executed
            if below_prot_mem1.trans != -1:
                if x_check < trap_locations[0] or x_check > trap_locations[1]:
                    return None
                elif y_check < trap_locations[0] or y_check > trap_locations[1]:
                    return None

        index += 1

    # if a complex passed all the previous checks we can turn and return a guide for the 'turn' function
    if prot.direction in ['N', 'S']:
        if rnd == 0:
            return turn_options[prot.direction][0], 'y', ['+', '-']
        else:
            return turn_options[prot.direction][1], 'y', ['-', '+']
    elif prot.direction in ['NE', 'SW']:
        if rnd == 0:
            return turn_options[prot.direction][0], 'y', ['-', '+']
        else:
            return turn_options[prot.direction][1], 'x', ['+', '-']
    elif prot.direction in ['E', 'W']:
        if rnd == 0:
            return turn_options[prot.direction][0], 'x', ['+', '-']
        else:
            return turn_options[prot.direction][1], 'x', ['-', '+']
    elif prot.direction in ['SE', 'NW']:
        if rnd == 0:
            return turn_options[prot.direction][0], 'y', ['+', '-']
        else:
            return turn_options[prot.direction][1], 'x', ['+', '-']


# this function perform the move of one protein for the turn function #
def move_to_turn(grid, protx, ax, distance, old_dir, new_dir):
    grid_len = len(grid)

    # from the moving guide perform a move on the x axis or y axis of the grid
    if ax == 'x':
        # remove the protein from the grid previous location
        grid[protx.location[0], protx.location[1], protx.location[2]] = 0
        # change the protein location in the protein location field
        protx.location = [adjust_xy(protx.location[0] + distance, grid_len), protx.location[1], protx.location[2]]
        # re-insert the protein to the grid in the new location
        grid[protx.location[0], protx.location[1], protx.location[2]] = protx
        # change the protein direction field
        if protx.direction not in [old_dir, inverse_direction(old_dir)]:
            raise SystemError("complex have more then two direction is found")
        if protx.direction == old_dir:
            protx.direction = new_dir
        else:
            if protx.direction == inverse_direction(old_dir):
                protx.direction = inverse_direction(new_dir)

    if ax == 'y':
        grid[protx.location[0], protx.location[1], protx.location[2]] = 0
        protx.location = [protx.location[0], adjust_xy(protx.location[1] + distance, grid_len), protx.location[2]]
        grid[protx.location[0], protx.location[1], protx.location[2]] = protx
        if protx.direction not in [old_dir, inverse_direction(old_dir)]:
            raise SystemError("complex have more then two direction is found")
        if protx.direction == old_dir:
            protx.direction = new_dir
        else:
            if protx.direction == inverse_direction(old_dir):
                protx.direction = inverse_direction(new_dir)

    return grid


# the main turn function - this is the only function that the 'Initiator.py' excess #
# this function choose a random protein and turn it's whole complex direction - if possible #
def turn(grid, all_protein, trap_locations):
    prot = random.choice(all_protein)  # chose a random protein
    if len(prot.complex) >= (2 * (len(grid) - 1)):  # if the zipper is very long (size of the grid's length) he won't
        # turn
        return grid

    # after choosing a protein, if he is a part of a complex, turn the complex from the middle.
    # this section change the chosen protein to the middle one in the complex
    if len(prot.complex) > 2:
        rnd = random.random()  # reduce the chance of a complex to turn by its length
        if rnd > (1 / (np.exp(len(prot.complex) - 6) + 1)):
            return grid
        if prot.direction in ['N', 'S']:
            lst, idx = sort_prot_by_direction(prot.complex.proteins, 'N')
        else:
            lst, idx = sort_prot_by_direction(prot.complex.proteins, 'W')
        middle_prot_idx = ((len(prot.complex) // 2) + idx) % len(prot.complex)
        if len(prot.complex) % 2 == 0:
            middle_prot_idx = random.choice([middle_prot_idx, middle_prot_idx - 1])
        prot = lst[middle_prot_idx]

    move_manual = can_complex_turn(grid, prot, trap_locations)  # check if the complex can turn
    if move_manual is None:
        return grid

    grid_len = len(grid)
    old_dir = prot.direction  # saves the old protein direction
    new_dir = move_manual[0]
    prot.complex.direction = [new_dir, inverse_direction(new_dir)]  # change the complex direction

    if old_dir in ['N', 'S']:
        index = 1
        while True:
            if grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 0] == 0 or \
                    grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 0] not in \
                    prot.complex.proteins:
                if grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 1] == 0 or \
                        grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 1] not in \
                        prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len), prot.location[1], 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1
        index = 1
        while True:
            if grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 0] == 0 or \
                    grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 0] not in \
                    prot.complex.proteins:
                if grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 1] == 0 or \
                        grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 1] not in \
                        prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len), prot.location[1], 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1

    if old_dir in ['W', 'E']:
        index = 1
        while True:
            if grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 0] == 0 or \
                    grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 0] not in \
                    prot.complex.proteins:
                if grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 1] == 0 or \
                        grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 1] not in \
                        prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 0] in \
                    prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[
                prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 1] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] + index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1
        index = 1
        while True:
            if grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 0] == 0 or \
                    grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 0] not in \
                    prot.complex.proteins:
                if grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 1] == 0 or \
                        grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 1] not in \
                        prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[
                prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 1] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[prot.location[0], adjust_xy(prot.location[1] - index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1

    if old_dir in ['NE', 'SW']:
        index = 1
        while True:
            if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] + index, grid_len),
                    0] == 0 or grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] + index,
                                                                                             grid_len), 0] not in \
                    prot.complex.proteins:
                if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] + index, grid_len),
                        1] == 0 or grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] +
                                                                                                 index, grid_len), 1] \
                        not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] + index, grid_len), 0] in \
                    prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 0], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 0], move_manual[1],
                                        -index, old_dir, new_dir)
            if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] + index, grid_len), 1] in \
                    prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 1], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 1], move_manual[1],
                                        -index, old_dir, new_dir)
            index += 1
        index = 1
        while True:
            if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] - index, grid_len),
                    0] == 0 or grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] - index,
                                                                                             grid_len), 0] not in \
                    prot.complex.proteins:
                if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] - index, grid_len),
                        1] == 0 or grid[adjust_xy(prot.location[0] + index, grid_len),
                                        adjust_xy(prot.location[1] - index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] - index, grid_len), 0] in \
                    prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 0], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 0], move_manual[1],
                                        -index, old_dir, new_dir)
            if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] - index, grid_len), 1] in \
                    prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 1], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 1], move_manual[1],
                                        -index, old_dir, new_dir)
            index += 1

    if old_dir in ['SE', 'NW']:
        index = 1
        while True:
            if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] - index, grid_len),
                    0] == 0 or grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] - index,
                                                                                             grid_len), 0] not in \
                    prot.complex.proteins:
                if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] - index, grid_len),
                        1] == 0 or grid[adjust_xy(prot.location[0] - index, grid_len),
                                        adjust_xy(prot.location[1] - index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] - index, grid_len),
                    0] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 0], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 0], move_manual[1],
                                        -index, old_dir, new_dir)
            if grid[adjust_xy(prot.location[0] - index, grid_len), adjust_xy(prot.location[1] - index, grid_len), 1] in \
                    prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 1], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] - index, grid_len),
                                                   adjust_xy(prot.location[1] - index, grid_len), 1], move_manual[1],
                                        -index, old_dir, new_dir)
                index += 1
        index = 1
        while True:
            if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] + index, grid_len),
                    0] == 0 or grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] + index,
                                                                                             grid_len), 0] not in \
                    prot.complex.proteins:
                if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] + index, grid_len),
                        1] == 0 or grid[adjust_xy(prot.location[0] + index, grid_len),
                                        adjust_xy(prot.location[1] + index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError("just a check for now, erase later")  # raise an error if the while loop had too
                # many iterations
            if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] + index, grid_len),
                    0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 0], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 0], move_manual[1],
                                        -index, old_dir, new_dir)
            if grid[adjust_xy(prot.location[0] + index, grid_len), adjust_xy(prot.location[1] + index, grid_len), 1] in \
                    prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 1], move_manual[1],
                                        index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[adjust_xy(prot.location[0] + index, grid_len),
                                                   adjust_xy(prot.location[1] + index, grid_len), 1], move_manual[1],
                                        -index, old_dir, new_dir)
            index += 1

    prot.direction = new_dir
    if prot.location[2] == 0:
        if grid[prot.location[0], prot.location[1], 1] in prot.complex.proteins:
            grid[prot.location[0], prot.location[1], 1].direction = inverse_direction(new_dir)
    if prot.location[2] == 1:
        if grid[prot.location[0], prot.location[1], 0] in prot.complex.proteins:
            grid[prot.location[0], prot.location[1], 0].direction = inverse_direction(new_dir)
    return grid


# create cis interaction  #


# the main create cis function - this is the only function that the 'Initiator.py' excess #
# this function choose a random protein and tries to create cis interaction based on it's direction #
def create_cis(grid, all_protein, all_complex, Cis_Mrx):
    grid_len = len(grid)
    prot = random.choice(all_protein)  # choose a random protein
    if prot.cis != -1:  # if the protein already have a cis interaction - don't do anything
        return grid, all_complex, all_protein
    if len(prot.complex) >= (2 * grid_len) - 2:  # to limit the max complex length to the grid length
        return grid, all_complex, all_protein

    dir = prot.direction
    x = prot.location[0]
    y = prot.location[1]
    direction_to_location = {'N': [1, 0], 'W': [0, 1], 'S': [-1, 0], 'E': [0, -1], 'SE': [-1, -1], 'SW': [-1, 1],
                             'NE': [1, -1], 'NW': [1, 1]}  # gives the coordination different from a cell to a nearby
    # cell in that direction calculate the new location
    loc_to_check = direction_to_location[dir]  # in order to find the cell the protein can interact with
    x_new = adjust_xy(x + loc_to_check[0], grid_len)
    y_new = adjust_xy(y + loc_to_check[1], grid_len)
    option_prot = grid[x_new, y_new, prot.location[2]]
    if option_prot == 0:  # if the optional cell is empty or already have other cis interaction - the action will not
        # be execute
        return grid, all_complex, all_protein
    if option_prot.cis != -1:
        return grid, all_complex, all_protein

    # in order to prevent intersection between complexes - this next 4 sections will check for it
    if prot.direction == 'SW':
        if grid[adjust_xy(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[adjust_xy(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].direction in \
                    ['NW', 'SE']:
                if len(grid[adjust_xy(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]].direction in \
                    ['NW', 'SE']:
                if len(grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if prot.direction == 'SE':
        if grid[adjust_xy(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[adjust_xy(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].direction in \
                    ['NE', 'SW']:
                if len(grid[adjust_xy(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]].direction in \
                    ['NE', 'SW']:
                if len(grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if prot.direction == 'NW':
        if grid[adjust_xy(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[adjust_xy(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].direction in \
                    ['SW', 'NE']:
                if len(grid[adjust_xy(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]].direction in \
                    ['SW', 'NE']:
                if len(grid[prot.location[0], adjust_xy(prot.location[1] + 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if prot.direction == 'NE':
        if grid[adjust_xy(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[adjust_xy(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].direction in \
                    ['SE', 'NW']:
                if len(grid[adjust_xy(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]].direction in \
                    ['SE', 'NW']:
                if len(grid[prot.location[0], adjust_xy(prot.location[1] - 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    # if the protein can interact in cis, we will use this chance function in order to maintain the inputted affinity
    if option_prot.direction == inverse_direction(prot.direction):
        # identify the cis affinity from the affinity matrix
        if not pd.isnull(Cis_Mrx.iloc[(prot.type - 1), (option_prot.type - 1)]):
            KA_cis = Cis_Mrx.iloc[(prot.type - 1), (option_prot.type - 1)]
        else:
            KA_cis = Cis_Mrx.iloc[(option_prot.type - 1), (prot.type - 1)]

        deltaG = np.exp(KA_cis) / (np.exp(KA_cis) + 1)  # apply the affinity function
        rnd = random.random()
        if rnd > deltaG:
            return grid, all_complex, all_protein

        # updating the complex and proteins based on the merge of both complexes
        new_prot_lst = option_prot.complex.proteins.copy()
        all_complex.remove(option_prot.complex)  # remove the old complex - since the two are merged into one
        for p in new_prot_lst:
            p.complex = prot.complex  # changing the complex field in each of the proteins to the new merged complex
        prot.complex.proteins.extend(new_prot_lst)  # merging the two protein lists in the merged complex

        # creating the new cis interaction
        prot.cis = option_prot
        option_prot.cis = prot

    return grid, all_complex, all_protein


# break cis interaction #

# this function return a list of all the proteins from the trans side of a specific protein #
def trans_side_tree(prot, lst):
    if prot.trans == -1:
        return lst
    lst.append(prot.trans)
    lst = cis_side_tree(prot.trans, lst)  # using recursion to keep on accessing the trans tree
    return lst


# this function return a list of all the proteins from the cis side of a specific protein #
def cis_side_tree(prot, lst):
    if prot.cis == -1:
        return lst
    lst.append(prot.cis)
    lst = trans_side_tree(prot.cis, lst)  # using recursion to keep on accessing the cis tree
    return lst


# the main break cis function - this is the only function that the 'Initiator.py' excess #
# this function choose a random protein and tries to break cis interaction if they have any #
def break_cis(grid, all_protein, all_complex, Cis_Mrx, frame, step):
    prot = random.choice(all_protein)  # choose a random protein
    if prot.cis == -1:  # if the protein don't have cis - it can't perform the break
        return grid, all_complex, all_protein

    # if the protein has cis interaction, we will use this chance function in order to maintain the inputted affinity
    # identify the cis affinity from the affinity matrix
    if not pd.isnull(Cis_Mrx.iloc[(prot.type - 1), (prot.cis.type - 1)]):
        KA_cis = Cis_Mrx.iloc[(prot.type - 1), (prot.cis.type - 1)]
    else:
        KA_cis = Cis_Mrx.iloc[(prot.cis.type - 1), (prot.type - 1)]

    deltaG = 1 / (np.exp(KA_cis) + 1)
    rnd = random.random()
    if rnd > (deltaG * (len(all_protein) / 2) / (len(grid) ** 2)):
        return grid, all_complex, all_protein

    new_prot_lst = cis_side_tree(prot, [])  # accessing the cis side of the interaction in order to start a new
    # complex after the break
    new_comp_name = "comp" + str(frame) + str(step)  # creating the new complex
    direction = prot.complex.direction
    locals()[new_comp_name] = Complex(new_prot_lst, direction)
    all_complex.append(locals()[new_comp_name])

    # removing the new complex proteins from the old complex protein list field
    for p in new_prot_lst:
        p.complex = locals()[new_comp_name]
        prot.complex.proteins.remove(p)

    # breaking the cis interaction from both proteins
    prot.cis.cis = -1
    prot.cis = -1

    return grid, all_complex, all_protein


# create trans interaction  #

# the main create trans function - this is the only function that the 'Initiator.py' excess #
# this function choose a random protein and tries to create trans interaction if they can occur #
def create_trans(grid, all_protein, all_complex, KAtransSelf, KAtransOther, trap_locations):
    prot = random.choice(all_protein)  # choose a random protein
    if len(prot.complex) >= (2 * len(grid) - 2):  # to limit the max complex length
        return grid, all_complex

    # check if we are in the diffusion trap area, if not we can't preform trans interactions
    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
        return grid, all_complex
    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
        return grid, all_complex
    # check if the protein already has a trans interaction
    if prot.trans != -1:
        return grid, all_complex

    # locating the optional protein for trans
    if prot.location[2] == 0:
        opt_prot = grid[prot.location[0], prot.location[1], 1]
    elif prot.location[2] == 1:
        opt_prot = grid[prot.location[0], prot.location[1], 0]
    if opt_prot == 0:  # if no optional protein is present, don't execute
        return grid, all_complex
    if opt_prot.trans != -1:
        raise SystemError("too many isoforms in one location")

    # if the protein don't have trans interaction, we will use this chance function in order to maintain the inputted
    # affinity
    if opt_prot.direction == inverse_direction(prot.direction):  # check if the direction of both proteins is good
        # for trans
        if prot.type == opt_prot.type:
            deltaG = np.exp(KAtransSelf) / (np.exp(KAtransSelf) + 1)
        else:
            deltaG = np.exp(KAtransOther) / (np.exp(KAtransOther) + 1)
        rnd = random.random()
        if rnd > deltaG:
            return grid, all_complex

        # merging the complexes
        new_prot_lst = opt_prot.complex.proteins.copy()
        all_complex.remove(opt_prot.complex)  # removing one of the old complexes
        for p in new_prot_lst:
            p.complex = prot.complex  # updating the proteins complex to the new merged one
        prot.complex.proteins.extend(new_prot_lst)  # updating the proteins list in the merged complex

        # creating the trans interaction in both proteins
        prot.trans = opt_prot
        opt_prot.trans = prot

    return grid, all_complex


# break trans interaction #

# the main break trans function - this is the only function that the 'Initiator.py' excess #
# this function choose a random protein and tries to break trans interaction if they present #
def break_trans(grid, all_protein, all_complex, KDtransSelf, KDtransOther, frame, step):
    prot = random.choice(all_protein)  # choose a random protein
    if prot.trans == -1:  # if the protein doesn't have trans interaction it will not continue
        return grid, all_complex

    # if the protein have trans interaction, we will use this chance function in order to maintain the inputted affinity
    if prot.type == prot.trans.type:
        deltaG = 1 / (np.exp(KDtransSelf) + 1)
    else:
        deltaG = 1 / (np.exp(KDtransOther) + 1)
    rnd = random.random()
    if rnd > (deltaG * (len(all_protein) / 2) / (len(grid) ** 2)):
        return grid, all_complex

    new_prot_lst = trans_side_tree(prot, [])  # accessing the trans side of the interaction in order to start a new
    # complex after the break
    new_comp_name = "comp" + str(frame) + str(step)  # creating the new complex
    direction = prot.complex.direction
    locals()[new_comp_name] = Complex(new_prot_lst, direction)
    all_complex.append(locals()[new_comp_name])

    # removing the new complex proteins from the old complex protein list field
    for p in new_prot_lst:
        p.complex = locals()[new_comp_name]
        prot.complex.proteins.remove(p)

    # breaking the trans interaction from both proteins
    prot.trans.trans = -1
    prot.trans = -1

    return grid, all_complex


########################################################################################################################

# functions for data extraction #


# saving the grid during experiment running for the 'vis_complex_fromfile' later (from Analysis.py) #
def save_grid(all_protein, all_complex, exp, frame, name):
    save_grid_file = open(name + "/Exp" + str(exp) + "/Raw_data/grid_frame" + str(frame), "x")
    save_grid_file.write("# file for visualization of the frame #\n" +
                         "# x location, y location, colorID, isTrans, complexID, isCis, membrane, isoform\n")  # adding
    # 2 lines of header

    prot_count = 0  # to make sure all proteins has been saved
    id_count = max([i.col for i in all_complex]) + 1  # find the current max ID and start assigning new ones from there

    for cx in all_complex:
        if len(cx) > 2:
            if cx.col == 0:
                cx.col = id_count
                id_count += 1
            for pr in cx.proteins:
                if pr.trans != -1:
                    if pr.cis != -1:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",True," + str(cx.col) + ",True," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",True," + str(cx.col) + ",False," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                else:
                    if pr.cis != -1:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",False," + str(cx.col) + ",True," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",False," + str(cx.col) + ",False," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                prot_count += 1
        else:
            if cx.col != 0:
                cx.col = 0
            for pr in cx.proteins:
                if pr.trans != -1:
                    if pr.trans.type != pr.type:  # check that the trans dimer consist of the same isoforms
                        print("trans dimer of protein with different isoforms has been identify. frame:", frame,
                              "exp:", exp)
                        print(pr.trans.type, pr.type)
                    if pr.cis != -1:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",True," + str(cx.col) + ",True," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",True," + str(cx.col) + ",False," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                else:
                    if pr.cis != -1:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",False," + str(cx.col) + ",True," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",False," + str(cx.col) + ",False," + str(pr.location[2]) + "," +
                                             str(pr.type) + "\n")
            prot_count += 1
    save_grid_file.close()


# analyse the grid basic stats for the main data files #
def analyse_grid(file_path1, file_path2, file_path3, all_protein, all_complex, grid):
    # counting the number of cis and trans interactions * 2, will be divided later when writing the files
    cis_amount = 0
    trans_amount = 0
    for p in all_protein:
        if p.cis != -1:
            cis_amount += 1
        if p.trans != -1:
            trans_amount += 1

    # creating a list with all the complexes sizes
    complex_size = list()
    for c in all_complex:
        length = len(c)
        complex_size.append(length)

    # creating a list with all the complexes larger then 3 proteins
    comp_size_3 = list()
    for i in complex_size:
        if i > 2:
            comp_size_3.append(i)

    # calculating the average size of complexes larger then 3
    if comp_size_3 == []:
        comp_avg = 0
    else:
        comp_avg = sum(comp_size_3) / len(comp_size_3)
    # finding the maximum complex length (even if smaller then 3)
    comp_max = max(complex_size)

    data = None
    try:
        data = open(file_path1, 'a')  # opening for append the 'Stats' csv file
        # writing the number of cis, trans, avg complex size and max complex size
        data.write(str(cis_amount / 2) + "," + str(trans_amount / 2) + "," + str(comp_avg) + "," + str(comp_max) + "\n")
    finally:
        if data is not None:
            data.close()

    full_comp_data = None
    try:
        full_comp_data = open(file_path2, 'a')  # opening for append the 'Full_Complex_Length' csv file
        for i in complex_size:  # writing all complexes length (even smaller than 3)
            full_comp_data.write(str(i) + ",")
        full_comp_data.write("\n")
    finally:
        if full_comp_data is not None:
            full_comp_data.close()

    # analyse the clusters of trans interactions. the result will be a list of clusters, each cluster will hold the
    # cell locations of proteins belong to the cluster. the analysis will be performed on membrane 0 since we only
    # cluster trans interactions, therefore the result should be the same#
    grid_len = len(grid)
    all_clusters = list()
    for i in range(grid_len):
        for j in range(grid_len):
            if grid[i, j, 0] != 0 and grid[i, j, 0].trans != -1 and grid[i, j, 0].flag == False:  # make sure we won't
                # check the same protein more than once
                cur_cluster = cluster(i, j, grid)
                if len(cur_cluster) > 1:
                    all_clusters.append(cur_cluster)

    full_cluster_data = None
    try:
        full_cluster_data = open(file_path3, 'a')  # opening for append the 'Full_Complex_Length' csv file
        for cl in all_clusters:
            for p in cl:
                # in the file each line is a frame. in each line clusters are separated by "," and each protein in a
                # cluster is separated by ";" and the x&y values are separated by ":"
                full_cluster_data.write(str(p[0]) + ":" + str(p[1]) + ";")
            full_cluster_data.write(",")
        full_cluster_data.write("\n")
    finally:
        if full_cluster_data is not None:
            full_cluster_data.close()

    for p in all_protein:  # after analysis all proteins get a correct, unchecked flag until the next frame
        p.flag = False
