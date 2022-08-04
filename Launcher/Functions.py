import numpy as np
import random
import math
import sys
import os
from zipfile import ZipFile


# class info #

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
        self.direction = directions
        self.proteins = prot_lst
        self.col = 0  # for visualization through save_grid function

    def __len__(self):
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
        self.location = location
        self.type = iso_type
        self.cis = cis
        self.trans = -1
        self.direction = direction
        self.complex = complex
        if location[2] == 0:  # for visualization in the heatmap
            self.num = 1
        if location[2] == 1:
            self.num = 4

    def __repr__(self):  # name in the grid - can be modified for debug
        return self.direction

    def __add__(self, other):  # for visualization in the heatmap
        if isinstance(other, int):
            return self.num + other
        if isinstance(other, Protein):
            return self.num + other.num

    def __radd__(self, other):  # for visualization in the heatmap
        if isinstance(other, int):
            return self.num + other
        if isinstance(other, Protein):
            return self.num + other.num


# create a grid #

def create_empty_grid(grid_len):  # create a grid with zeros
    grid = np.zeros((grid_len, grid_len, 2), dtype=Protein)
    return grid


def taking_care_of_edges(x, y, grid_len):  # correct x,y in case of getting to the edge
    if x < 0:
        x = grid_len - 1
    if x > grid_len - 1:
        x = 0
    if y < 0:
        y = grid_len - 1
    if y > grid_len - 1:
        y = 0
    return x, y


def find_spot_for_dimer(grid, membrane):  # locate spots for dimer insertion while creating the initial grid
    grid_len = len(grid)
    x1 = random.choice(range(grid_len))
    y1 = random.choice(range(grid_len))
    while grid[x1, y1, membrane] != 0:
        x1 = random.choice(range(grid_len))
        y1 = random.choice(range(grid_len))
    direction = random.choice([['N', 'S'], ['W', 'E'], ['SE', 'NW'], ['SW', 'NE']])

    if direction == ['N', 'S']:
        x2 = x1 + 1
        y2 = y1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        x2 = x1 - 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)  # using recursion if the spot wasn't free

    if direction == ['W', 'E']:
        x2 = x1
        y2 = y1 + 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        y2 = y1 - 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)  # using recursion if the spot wasn't free

    if direction == ['SE', 'NW']:
        x2 = x1 - 1
        y2 = y1 - 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x_edge_correct(x1 - 1, grid_len), y1, membrane] != 0:  # make sure complexes won't intersect
                if grid[x_edge_correct(x1 - 1, grid_len), y1, membrane].direction in ['NE', 'SW']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        x2 = x1 + 1
        y2 = y1 + 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x_edge_correct(x1 + 1, grid_len), y1, membrane] != 0:  # make sure complexes won't intersect
                if grid[x_edge_correct(x1 + 1, grid_len), y1, membrane].direction in ['NE', 'SW']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)  # using recursion if the spot wasn't free

    if direction == ['SW', 'NE']:
        x2 = x1 - 1
        y2 = y1 + 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x_edge_correct(x1 - 1, grid_len), y1, membrane] != 0:  # make sure complexes won't intersect
                if grid[x_edge_correct(x1 - 1, grid_len), y1, membrane].direction in ['NW', 'SE']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x1, y1, membrane], [x2, y2, membrane]], direction
        x2 = x1 + 1
        y2 = y1 - 1
        x2, y2 = taking_care_of_edges(x2, y2, grid_len)
        if grid[x2, y2, membrane] == 0:
            if grid[x_edge_correct(x1 + 1, grid_len), y1, membrane] != 0:  # make sure complexes won't intersect
                if grid[x_edge_correct(x1 + 1, grid_len), y1, membrane].direction in ['NW', 'SE']:
                    return find_spot_for_dimer(grid, membrane)
            return [[x2, y2, membrane], [x1, y1, membrane]], direction
        else:
            return find_spot_for_dimer(grid, membrane)  # using recursion if the spot wasn't free

    raise SystemError("dimer skipped from insertion to grid")


def insert_protein_to_grid(grid, prot):
    grid[prot.location[0], prot.location[1], prot.location[2]] = prot
    return grid


def init_grid(grid_len, isoAmount):  # create a grid with dimers inserted randomly
    grid = create_empty_grid(grid_len)

    membrane0_iso_lst = [1]
    membrane1_iso_lst = [1]
    membrane0_iso_dict = dict()  # creating isoforms inventory for choosing later
    for i in membrane0_iso_lst:
        membrane0_iso_dict[i] = isoAmount
    membrane1_iso_dict = dict()
    for i in membrane1_iso_lst:
        membrane1_iso_dict[i] = isoAmount

    all_protein_id = range(1, (isoAmount * 2) + 1)
    all_protein = ['prot' + str(i) for i in all_protein_id]  # list of all the proteins
    all_complex_id = range(1, (isoAmount) + 1)
    all_complex = ['comp' + str(i) for i in all_complex_id]  # list of all the complexes

    prot_index = 0
    for comp_index in range((len(all_complex)) // 2):
        comp_name = all_complex[comp_index]
        loc = find_spot_for_dimer(grid, 0)
        prot1_name = all_protein[prot_index]
        prot2_name = all_protein[prot_index + 1]
        typ1 = random.choice(membrane0_iso_lst)
        while membrane0_iso_dict[typ1] == 0:
            typ1 = random.choice(membrane0_iso_lst)
        membrane0_iso_dict[typ1] -= 1
        typ2 = random.choice(membrane0_iso_lst)
        while membrane0_iso_dict[typ2] == 0:
            typ2 = random.choice(membrane0_iso_lst)
        membrane0_iso_dict[typ2] -= 1

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

    for comp_index in range((len(all_complex)) // 2, len(all_complex)):
        comp_name = all_complex[comp_index]
        loc = find_spot_for_dimer(grid, 1)
        prot1_name = all_protein[prot_index]
        prot2_name = all_protein[prot_index + 1]
        typ1 = random.choice(membrane1_iso_lst)
        while membrane1_iso_dict[typ1] == 0:
            typ1 = random.choice(membrane1_iso_lst)
        membrane1_iso_dict[typ1] -= 1
        typ2 = random.choice(membrane1_iso_lst)
        while membrane1_iso_dict[typ2] == 0:
            typ2 = random.choice(membrane1_iso_lst)
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
    return grid, all_protein, all_complex


# move #

def sort_prot_by_direction(prot_lst,
                           direction):  # sorting the protein list by the move direction. from first to move to last to move
    index = 0
    if direction in ['N', 'NW', 'NE']:
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[0])
        for i in range(len(sorted_lst) - 1):  # change the start index in case of complexes that are crossing the edges
            if sorted_lst[i + 1].location[0] - sorted_lst[i].location[0] not in [0, 1]:
                index = i + 1
    elif direction in ['S', 'SW', 'SE']:
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[0], reverse=True)
        for i in range(len(sorted_lst) - 1):  # change the start index in case of complexes that are crossing the edges
            if sorted_lst[i].location[0] - sorted_lst[i + 1].location[0] not in [0, 1]:
                index = i + 1
    elif direction == 'E':
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[1], reverse=True)
        for i in range(len(sorted_lst) - 1):  # change the start index in case of complexes that are crossing the edges
            if sorted_lst[i].location[1] - sorted_lst[i + 1].location[1] not in [0, 1]:
                index = i + 1
    elif direction == 'W':
        sorted_lst = sorted(prot_lst, key=lambda prot: prot.location[1])
        for i in range(len(sorted_lst) - 1):  # change the start index in case of complexes that are crossing the edges
            if sorted_lst[i + 1].location[1] - sorted_lst[i].location[1] not in [0, 1]:
                index = i + 1
    return sorted_lst, index


def check_one_move(grid, prot, direction, same_dir=False):  # check if one protein can move in a direction. same_dir
    # should be 'True' if we move to a direction that is in the direction of the complex - to prevent complex
    # intersecting
    grid_len = len(grid)
    side_to_loc = {'N': [-1, 0], 'S': [1, 0], 'W': [0, -1], 'E': [0, 1], 'NE': [-1, 1], 'NW': [-1, -1], 'SE': [1, 1],
                   'SW': [1, -1]}
    dif_x = side_to_loc[direction][0]
    dif_y = side_to_loc[direction][1]
    x_check = prot.location[0] + dif_x
    y_check = prot.location[1] + dif_y
    x_check, y_check = taking_care_of_edges(x_check, y_check, grid_len)
    if grid[x_check, y_check, prot.location[2]] == 0:
        if not same_dir:
            if prot.direction in ['N', 'E', 'S', 'W']:
                return True

            else:
                if prot.location[2] == 0:
                    if grid[prot.location[0], prot.location[1], 1] == 0:
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

            if direction == 'NW':
                if grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]] == 0:
                    if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]] == 0:
                        return True

            if direction == 'NE':
                if grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]] == 0:
                    if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]] == 0:
                        return True

            if direction == 'SW':
                if grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]] == 0:
                    if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]] == 0:
                        return True

            if direction == 'SE':
                if grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]] == 0:
                    if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]] == 0:
                        return True
    else:
        return False


def move_one(grid, prot, direction):  # move one protein in a direction
    grid_len = len(grid)
    side_to_loc = {'N': [-1, 0], 'S': [1, 0], 'W': [0, -1], 'E': [0, 1], 'NE': [-1, 1], 'NW': [-1, -1], 'SE': [1, 1],
                   'SW': [1, -1]}
    x_new = prot.location[0] + side_to_loc[direction][0]
    y_new = prot.location[1] + side_to_loc[direction][1]
    x_new, y_new = taking_care_of_edges(x_new, y_new, grid_len)
    z_new = prot.location[2]
    old_loc = prot.location
    prot.location = [x_new, y_new, z_new]
    grid[x_new, y_new, z_new] = prot
    grid[old_loc[0], old_loc[1], old_loc[2]] = 0
    return grid


def move(grid, all_protein, trap_locations):  # move complex to a random location, if possible
    ran_prot = random.choice(all_protein)
    complex = ran_prot.complex
    if len(complex) >= (2*(len(grid) - 1)):  # if the zipper is very long he won't move
        return grid
    if not isinstance(complex, Complex):
        raise SystemError("trying to move a non complex type")
    direction = random.choice(['N', 'S', 'W', 'E', 'NW', 'NE', 'SW', 'SE'])

    side_to_loc = {'N': [-1, 0], 'S': [1, 0], 'W': [0, -1], 'E': [0, 1], 'NE': [-1, 1], 'NW': [-1, -1], 'SE': [1, 1],
                   'SW': [1, -1]}
    for p in complex.proteins:  # check that a complex with trans does not cross the diffusion trap
        if p.trans != -1:
            x_new = p.location[0] + side_to_loc[direction][0]
            y_new = p.location[1] + side_to_loc[direction][1]
            x_new, y_new = taking_care_of_edges(x_new, y_new, len(grid))
            if x_new < trap_locations[0] or x_new > trap_locations[1]:
                return grid
            if y_new < trap_locations[0] or y_new > trap_locations[1]:
                return grid

    if direction in complex.direction:
        sorted_prot_lst, index = sort_prot_by_direction(complex.proteins, direction)
        if check_one_move(grid, sorted_prot_lst[index], direction, True):
            if sorted_prot_lst[index].trans != -1:  # check the trans of the first protein in the chain
                if check_one_move(grid, sorted_prot_lst[index].trans, direction, True):
                    for i in range(index, len(sorted_prot_lst)):
                        grid = move_one(grid, sorted_prot_lst[i], direction)
                    for i in range(0, index):
                        grid = move_one(grid, sorted_prot_lst[i], direction)
            if sorted_prot_lst[index].trans == -1:
                if sorted_prot_lst[index].cis != -1:  # check the trans of the cis to the first protein in the chain
                    if sorted_prot_lst[index].cis.trans != -1:
                        if check_one_move(grid, sorted_prot_lst[index].cis.trans, direction, True):
                            for i in range(index, len(sorted_prot_lst)):
                                grid = move_one(grid, sorted_prot_lst[i], direction)
                            for i in range(0, index):
                                grid = move_one(grid, sorted_prot_lst[i], direction)

    else:
        for prot in complex.proteins:
            if not check_one_move(grid, prot, direction, False):
                return grid
        for prot in complex.proteins:
            grid = move_one(grid, prot, direction)
    return grid


# turn #

def x_edge_correct(x, grid_len):  # correct the x values of the edges
    if x < 0:
        x = grid_len + x
    if x >= grid_len:
        x = x - grid_len
    return x


def y_edge_correct(y, grid_len):  # correct the y values of the edges
    if y < 0:
        y = grid_len + y
    if y >= grid_len:
        y = y - grid_len
    return y


def can_complex_turn(grid, prot, trap_locations):  # check if a complex can more around a chosen protein (right or left - randomly chosen)
    turn_options = {'N': ['NE', 'NW'], 'S': ['SW', 'SE'], 'E': ['SE', 'NE'], 'W': ['NW', 'SW'], 'SE': ['S', 'E'],
                    'NE': ['N', 'E'], 'SW': ['S', 'W'], 'NW': ['N', 'W']}
    rnd = random.choice([0, 1])
    grid_len = len(grid)

    if prot.direction in ['N', 'S'] and rnd == 0:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + i, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] + index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] + index, grid_len) > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + i,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] + index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] + index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] == 0 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] == 0 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - i, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] - index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] - index, grid_len) > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - i,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] - index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] - index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][0], 'y', ['+', '-']

    if prot.direction in ['N', 'S'] and rnd == 1:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - i,
                                                                                               grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] - index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] - index, grid_len) > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - i,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index + 1, grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] - index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] - index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] == 0 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] == 0 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + i,
                                                                                               grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] + index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] + index, grid_len) > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + i,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index - 1, grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if y_edge_correct(prot.location[1] + index, grid_len) < trap_locations[0] or y_edge_correct(prot.location[1] + index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][1], 'y', ['-', '+']

    if prot.direction in ['NE', 'SW'] and rnd == 0:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index - i, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + 1,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index - i,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + 1,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index + i,
                                                                                           grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - 1,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index + i,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - 1,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][0], 'y', ['-', '+']

    if prot.direction in ['NE', 'SW'] and rnd == 1:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] - index + i, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                               grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] + index,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] - index + i, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] + index,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index - i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                               grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] - index,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index - i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] - index,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][1], 'x', ['+', '-']

    if prot.direction in ['E', 'W'] and rnd == 0:
        index = 1
        while True:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] not in prot.complex.proteins:
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + i, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] + index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] + index, grid_len) > trap_locations[1]:
                        return None
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + i, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] + index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] + index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] not in prot.complex.proteins:
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - i, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] - index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] - index, grid_len) > trap_locations[1]:
                        return None
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] - index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] - index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][0], 'x', ['+', '-']

    if prot.direction in ['E', 'W'] and rnd == 1:
        index = 1
        while True:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] not in prot.complex.proteins:
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - i, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] - index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] - index, grid_len) > trap_locations[1]:
                        return None
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - i, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index + 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] - index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] - index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] not in prot.complex.proteins:
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] == 0 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 0] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] + index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] + index, grid_len) > trap_locations[1]:
                        return None
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index - 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].complex) > 2:
                            return None
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if x_edge_correct(prot.location[0] + index, grid_len) < trap_locations[0] or x_edge_correct(prot.location[0] + index, grid_len) > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][1], 'x', ['-', '+']

    if prot.direction in ['SE', 'NW'] and rnd == 0:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1] == 0 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                print(index, prot, prot.location)
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index + i, grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - 1,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index + i,
                                                                                           grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - 1,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index - i,
                                                                                           grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + 1, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + 1,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index - i, grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + 1, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + 1,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][0], 'y', ['+', '-']

    if prot.direction in ['SE', 'NW'] and rnd == 1:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] - index + i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                               grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] - index,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] - index + i, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] - 1, grid_len), y_edge_correct(prot.location[1] - index,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index - i, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                               grid_len), 0] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] + index,grid_len), 1].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                for i in range(1, index + 1):
                    if grid[
                        x_edge_correct(prot.location[0] + index - i, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                               grid_len), 1] != 0:
                        return None
                if index > 1 or grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans == -1:
                    if grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 0] != 0:  # check the other membrane to prevent intersections between complexes
                        if len(grid[x_edge_correct(prot.location[0] + 1, grid_len), y_edge_correct(prot.location[1] + index,grid_len), 0].complex) > 2:
                            return None
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index, grid_len), 1].trans != -1:  # check that a trans protein won't leave the trap while turning
                    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:
                        return None
            index += 1
        return turn_options[prot.direction][1], 'x', ['+', '-']


def inverse_direction(direction):  # find the inverse direction for a given direction (180)
    dic = {'N': 'S', 'S': 'N', 'E': 'W', 'W': 'E', 'NW': 'SE', 'NE': 'SW', 'SE': 'NW', 'SW': 'NE'}
    inve_dir = dic[direction]
    return inve_dir


def move_to_turn(grid, protx, ax, distance, old_dir, new_dir):
    grid_len = len(grid)
    if ax == 'x':
        grid[protx.location[0], protx.location[1], protx.location[2]] = 0
        protx.location = [x_edge_correct(protx.location[0] + distance, grid_len), protx.location[1], protx.location[2]]
        grid[protx.location[0], protx.location[1], protx.location[2]] = protx
        if protx.direction not in [old_dir, inverse_direction(old_dir)]:
            raise SystemError("complex have more then two direction is found")
        if protx.direction == old_dir:
            protx.direction = new_dir
        else:
            if protx.direction == inverse_direction(old_dir):
                protx.direction = inverse_direction(new_dir)
    if ax == 'y':
        grid[protx.location[0], protx.location[1], protx.location[2]] = 0
        protx.location = [protx.location[0], y_edge_correct(protx.location[1] + distance, grid_len), protx.location[2]]
        grid[protx.location[0], protx.location[1], protx.location[2]] = protx
        if protx.direction not in [old_dir, inverse_direction(old_dir)]:
            raise SystemError("complex have more then two direction is found")
        if protx.direction == old_dir:
            protx.direction = new_dir
        else:
            if protx.direction == inverse_direction(old_dir):
                protx.direction = inverse_direction(new_dir)
    return grid


def turn(grid, all_protein, trap_locations):  # run the check functions and turn the complex if possible
    prot = random.choice(all_protein)
    if len(prot.complex) >= (2*(len(grid) - 1)):  # if the zipper is very long he won't move
        return grid
    ## test for turning from the middle of the complex only ##
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

    move_manual = can_complex_turn(grid, prot, trap_locations)
    if move_manual is None:
        return grid
    grid_len = len(grid)
    old_dir = prot.direction
    new_dir = move_manual[0]
    prot.complex.direction = [new_dir, inverse_direction(new_dir)]

    if old_dir in ['N', 'S']:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] == 0 or grid[
                x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] == 0 or grid[
                    x_edge_correct(prot.location[0] - index, grid_len), prot.location[
                        1], 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), prot.location[1], 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), prot.location[
                        1], 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[
                x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[
                        x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[
                        x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[
                x_edge_correct(prot.location[0] + index, grid_len), prot.location[1], 1] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[
                                            1], 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] + index, grid_len), prot.location[
                                            1], 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1

    if old_dir in ['W', 'E']:
        index = 1
        while True:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] == 0 or grid[
                prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] not in prot.complex.proteins:
                if grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] == 0 or grid[
                    prot.location[0], y_edge_correct(prot.location[1] + index,
                                                     grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[
                prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[
                prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[prot.location[0], y_edge_correct(prot.location[1] + index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1
        index = 1
        while True:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] == 0 or grid[
                prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 0] not in prot.complex.proteins:
                if grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] == 0 or grid[
                    prot.location[0], y_edge_correct(prot.location[1] - index,
                                                     grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[prot.location[0], y_edge_correct(prot.location[1] - index,
                                                     grid_len), 0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[prot.location[0], y_edge_correct(prot.location[1] - index,
                                                                                    grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[prot.location[0], y_edge_correct(prot.location[1] - index,
                                                                                    grid_len), 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[
                prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid,
                                        grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid,
                                        grid[prot.location[0], y_edge_correct(prot.location[1] - index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1

    if old_dir in ['NE', 'SW']:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] + index, grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] + index, grid_len), 0], move_manual[1], -index, old_dir,
                                        new_dir)
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] + index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] + index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(
                                            prot.location[1] - index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(
                                            prot.location[1] - index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
            index += 1

    if old_dir in ['SE', 'NW']:
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] - index, grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] - index, grid_len), 0], move_manual[1], -index, old_dir,
                                        new_dir)
            if grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(prot.location[1] - index,
                                                                                       grid_len), 1] in prot.complex.proteins:
                if move_manual[2][0] == '+':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] - index, grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][0] == '-':
                    grid = move_to_turn(grid,
                                        grid[x_edge_correct(prot.location[0] - index, grid_len), y_edge_correct(
                                            prot.location[1] - index, grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
                index += 1
        index = 1
        while True:
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] == 0 or grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                   grid_len), 0] not in prot.complex.proteins:
                if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 1] == 0 or grid[
                    x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 1] not in prot.complex.proteins:
                    break
            if index > grid_len:
                raise SystemError(
                    "just a check for now, erase later")  # just to check everything is okay with the loop...
                # return None  # if the complex is in the length of the grid, don't turn
            if grid[x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                       grid_len), 0] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid, grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 0],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid, grid[
                        x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                           grid_len), 0],
                                        move_manual[1], -index, old_dir, new_dir)
            if grid[
                x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(prot.location[1] + index,
                                                                                   grid_len), 1] in prot.complex.proteins:
                if move_manual[2][1] == '+':
                    grid = move_to_turn(grid,
                                        grid[
                                            x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(
                                                prot.location[1] + index,
                                                grid_len), 1],
                                        move_manual[1], index, old_dir, new_dir)
                if move_manual[2][1] == '-':
                    grid = move_to_turn(grid,
                                        grid[
                                            x_edge_correct(prot.location[0] + index, grid_len), y_edge_correct(
                                                prot.location[1] + index,
                                                grid_len), 1],
                                        move_manual[1], -index, old_dir, new_dir)
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

def create_cis(grid, all_protein, all_complex, KAcisSelf, KAcisOther):  # create cis interactions between two proteins
    grid_len = len(grid)
    prot = random.choice(all_protein)
    if prot.cis != -1:
        return grid, all_complex, all_protein
    if len(prot.complex) >= (2 * grid_len) - 2:  # to limit the max complex length
        return grid, all_complex, all_protein
    dir = prot.direction
    x = prot.location[0]
    y = prot.location[1]
    direction_to_location = {'N': [1, 0], 'W': [0, 1], 'S': [-1, 0], 'E': [0, -1], 'SE': [-1, -1], 'SW': [-1, 1],
                             'NE': [1, -1], 'NW': [1, 1]}
    loc_to_check = direction_to_location[dir]
    x_new = x_edge_correct(x + loc_to_check[0], grid_len)
    y_new = y_edge_correct(y + loc_to_check[1], grid_len)
    option_prot = grid[x_new, y_new, prot.location[2]]
    if option_prot == 0:
        return grid, all_complex, all_protein
    if option_prot.cis != -1:
        return grid, all_complex, all_protein

    if prot.direction == 'SW':  # in order to prevent intersections between complexes
        if grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].direction in ['NW', 'SE']:
                if len(grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]].direction in ['NW', 'SE']:
                if len(grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if prot.direction == 'SE':  # in order to prevent intersections between complexes
        if grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
                 return grid, all_complex, all_protein
            if grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].direction in ['NE', 'SW']:
                if len(grid[x_edge_correct(prot.location[0] - 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]].direction in ['NE', 'SW']:
                if len(grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if prot.direction == 'NW':  # in order to prevent intersections between complexes
        if grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].direction in ['SW','NE']:
                if len(grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]].direction in ['SW','NE']:
                if len(grid[prot.location[0], y_edge_correct(prot.location[1] + 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if prot.direction == 'NE':  # in order to prevent intersections between complexes
        if grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
                return grid, all_complex, all_protein
            if grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].direction in ['SE', 'NW']:
                if len(grid[x_edge_correct(prot.location[0] + 1, grid_len), prot.location[1], prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein
        if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]] != 0:
            if grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]].direction in ['SE', 'NW']:
                if len(grid[prot.location[0], y_edge_correct(prot.location[1] - 1, grid_len), prot.location[2]].complex) > 2:
                    return grid, all_complex, all_protein

    if option_prot.direction == inverse_direction(prot.direction):
        if prot.type == option_prot.type:
            deltaG = np.exp(KAcisSelf) / (np.exp(KAcisSelf) + 1)
        else:
            deltaG = np.exp(KAcisOther) / (np.exp(KAcisOther) + 1)
        rnd = random.random()
        if rnd > deltaG:
            return grid, all_complex, all_protein

        new_prot_lst = option_prot.complex.proteins.copy()
        all_complex.remove(option_prot.complex)
        for p in new_prot_lst:
            p.complex = prot.complex
        prot.complex.proteins.extend(new_prot_lst)

        prot.cis = option_prot
        option_prot.cis = prot

    return grid, all_complex, all_protein


# break cis interaction #

def trans_side_tree(prot, lst):  # locate all the proteins on the trans side of a protein
    if prot.trans == -1:
        return lst
    lst.append(prot.trans)
    lst = cis_side_tree(prot.trans, lst)
    return lst


def cis_side_tree(prot, lst):  # locate all the proteins on the cis side of a protein
    if prot.cis == -1:
        return lst
    lst.append(prot.cis)
    lst = trans_side_tree(prot.cis, lst)
    return lst


def break_cis(grid, all_protein, all_complex, KDcisSelf, KDcisOther, frame,
              step):  # break a cis bond between two proteins
    prot = random.choice(all_protein)
    if prot.cis == -1:
        return grid, all_complex, all_protein
    if prot.type == prot.cis.type:
        deltaG = 1 / (np.exp(KDcisSelf) + 1)
    else:
        deltaG = 1 / (np.exp(KDcisOther) + 1)
    rnd = random.random()
    if rnd > (deltaG*(len(all_protein)/2)/(len(grid)**2)):
        return grid, all_complex, all_protein

    new_prot_lst = cis_side_tree(prot, [])
    new_comp_name = "comp" + str(frame) + str(step)
    direction = prot.complex.direction
    locals()[new_comp_name] = Complex(new_prot_lst, direction)
    all_complex.append(locals()[new_comp_name])

    for p in new_prot_lst:
        p.complex = locals()[new_comp_name]
        prot.complex.proteins.remove(p)

    prot.cis.cis = -1
    prot.cis = -1
    return grid, all_complex, all_protein


# create trans interaction  #

def create_trans(grid, all_protein, all_complex, KAtransSelf, KAtransOther, trap_locations):  # create a trans interaction between two proteins
    prot = random.choice(all_protein)
    if len(prot.complex) >= (2 * len(grid)) - 2:  # to limit the max complex length
        return grid, all_complex
    if prot.location[0] < trap_locations[0] or prot.location[0] > trap_locations[1]:  # check we are in the difussion trap
        return grid, all_complex
    if prot.location[1] < trap_locations[0] or prot.location[1] > trap_locations[1]:
        return grid, all_complex
    if prot.trans != -1:
        return grid, all_complex
    if prot.location[2] == 0:
        opt_prot = grid[prot.location[0], prot.location[1], 1]
    elif prot.location[2] == 1:
        opt_prot = grid[prot.location[0], prot.location[1], 0]
    if opt_prot == 0:
        return grid, all_complex
    if opt_prot.trans != -1:
        raise SystemError("too many isoforms in one location")
    if opt_prot.direction == inverse_direction(prot.direction):
        if prot.type == opt_prot.type:
            deltaG = np.exp(KAtransSelf) / (np.exp(KAtransSelf) + 1)
        else:
            deltaG = np.exp(KAtransOther) / (np.exp(KAtransOther) + 1)
        rnd = random.random()
        if rnd > deltaG:
            return grid, all_complex

        new_prot_lst = opt_prot.complex.proteins.copy()
        all_complex.remove(opt_prot.complex)
        for p in new_prot_lst:
            p.complex = prot.complex
        prot.complex.proteins.extend(new_prot_lst)

        prot.trans = opt_prot
        opt_prot.trans = prot

    return grid, all_complex


# break trans interaction #

def break_trans(grid, all_protein, all_complex, KDtransSelf, KDtransOther, frame,
                step):  # break trans interaction between two proteins
    prot = random.choice(all_protein)
    if prot.trans == -1:
        return grid, all_complex
    if prot.type == prot.trans.type:
        deltaG = 1 / (np.exp(KDtransSelf) + 1)
    else:
        deltaG = 1 / (np.exp(KDtransOther) + 1)
    rnd = random.random()
    if rnd > (deltaG*(len(all_protein)/2)/(len(grid)**2)):
        return grid, all_complex

    new_prot_lst = trans_side_tree(prot, [])
    new_comp_name = "comp" + str(frame) + str(step)
    direction = prot.complex.direction
    locals()[new_comp_name] = Complex(new_prot_lst, direction)
    all_complex.append(locals()[new_comp_name])

    for p in new_prot_lst:
        p.complex = locals()[new_comp_name]
        prot.complex.proteins.remove(p)

    prot.trans.trans = -1
    prot.trans = -1

    return grid, all_complex


# diffusion trap locations #

def dif_trap_loc(grid_len, dif_trap):  # create a list with [xy_min, xy_max] that reprasnt the start and end of the trap x,y
    trap_len = math.floor(math.sqrt(dif_trap * (grid_len ** 2)))
    xy_min = (grid_len - trap_len) / 2
    xy_max = grid_len - xy_min
    xy_min = math.floor(xy_min)
    xy_max = math.floor(xy_max)
    if not isinstance(xy_min, int) or not isinstance(xy_max, int):
        raise SystemError("could not find the trap locations")
    trap_locations = [xy_min, xy_max]
    return trap_locations


# non-uniform step function #
def StepNonUni(frame):
    if frame < 25:
        return 1000
    elif frame < 50:
        return 5000
    elif frame < 75:
        return 100000
    else:
        return 300000


# saving the grid during experiment running for vis_complex_fromfile later from Analysis.py#
def save_grid(all_protein, all_complex, exp, frame, name):
    save_grid_file = open(name + "/Exp" + str(exp) + "/Raw_data/grid_frame" + str(frame), "x")
    save_grid_file.write("# file for visualization of the frame #\n" +
                         "# x location, y location, colorID, isTrans, complexID, isCis, membrane\n")  # adding 2 lines of header

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
                                             ",True," + str(cx.col) + ",True," + str(pr.location[2]) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",True," + str(cx.col) + ",False," + str(pr.location[2]) + "\n")
                else:
                    if pr.cis != -1:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",False," + str(cx.col) + ",True," + str(pr.location[2]) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",c" + str(pr.type) +
                                             ",False," + str(cx.col) + ",False," + str(pr.location[2]) + "\n")
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
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",t" + str(pr.type) +
                                             ",True," + str(cx.col) + ",True," + str(pr.location[2]) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",t" + str(pr.type) +
                                             ",True," + str(cx.col) + ",False," + str(pr.location[2]) + "\n")
                else:
                    if pr.cis != -1:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",o" + str(pr.type) +
                                             ",False," + str(cx.col) + ",True," + str(pr.location[2]) + "\n")
                    else:
                        save_grid_file.write(str(pr.location[0]) + "," + str(pr.location[1]) + ",o" + str(pr.type) +
                                             ",False," + str(cx.col) + ",False," + str(pr.location[2]) + "\n")
            prot_count += 1
    save_grid_file.close()


# analyse a grid for data files #
def analyse_grid(file_path1, file_path2, all_protein, all_complex):
    cis_amount = 0
    trans_amount = 0
    for p in all_protein:
        if p.cis != -1:
            cis_amount += 1
        if p.trans != -1:
            trans_amount += 1

    complex_size = list()
    for c in all_complex:
        length = len(c)
        complex_size.append(length)

    comp_size_3 = list()
    for i in complex_size:
        if i > 2:
            comp_size_3.append(i)

    if comp_size_3 == []:
        comp_avg = 0
    else:
        comp_avg = sum(comp_size_3) / len(comp_size_3)
    comp_max = max(complex_size)

    data = None
    try:
        data = open(file_path1, 'a')
        data.write(str(cis_amount / 2) + "," + str(trans_amount / 2) + "," + str(comp_avg) + "," + str(comp_max) + "\n")
    finally:
        if data is not None:
            data.close()

    full_comp_data = None
    try:
        full_comp_data = open(file_path2, 'a')
        for i in complex_size:
            full_comp_data.write(str(i)+",")
        full_comp_data.write("\n")
    finally:
        if full_comp_data is not None:
            full_comp_data.close()
