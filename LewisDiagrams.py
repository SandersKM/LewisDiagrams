###############################
# CSCI 150 Fall, 2016
# Kate Sanders 
# Project 3 - Lewis Structures
###############################

import turtle

# self.atom_period_group_eneg is a dictionary containing
# all nonmetals with their respective period and group
# on the Periodic Table & electronegativity
# periods (rows) corresponds to the atom's orbital level
# groups (columns) corresponds to the number of
# valence electrons an atom has
# electronegativity is the amount a chemical property
# describing how well an atom can attract an electron to itself.
# These values are from:
# http://sciencenotes.org/list-of-electronegativity-values-of-the-elements/
# I gave H and Noble Gasses an eneg value of 100
# Noble Gasses don't have electronegativity numbers
# H's eneg value is normally low but it cannot be the central atom
atom_period_group_eneg = {"H": [1, 1, 100], "He": [1, 2, 100], \
                    "B": [2, 3, 2.04], "C": [2, 4, 2.55], \
                    "N": [2, 5, 3.04], "O": [2, 6, 3.44], \
                    "F": [2, 7, 3.98], "Ne": [2, 8, 100], \
                    "Si": [3, 4, 1.90], "P": [3, 5, 2.19], \
                    "S": [3, 6, 2.58], "Cl": [3, 7, 3.16], \
                    "Ar": [3, 8, 100], "As": [4, 5, 2.18], \
                    "Se": [4, 6, 2.55], "Br": [4, 7, 2.96],
                    "Kr": [4, 8, 3.00], "Te": [5, 6, 2.1], \
                    "I": [5, 7, 2.66], "Xe": [5, 8, 2.6], \
                    "At": [5, 9, 2.2], "Rn": [5, 10, 0.89]}

e_mol_geometry = {"linear": ["linear", "linear", "linear", "linear"], \
                  "trigonal planar": \
                  ["trigonal planar", "bent"], "tetrahedral": \
                  ["tetrahedral", "trigonal pyramidal", "bent"], \
                  "trigonal bipyramidal": \
                  ["trigonal bipyramidal", "seesaw", "t-shaped", \
                   "linear"], "octahedral": \
                   ["octahedral", "square pyramidal", "square planar"], \
                   "pentagonal bipyramidal": \
                   ["pentagonal bipyramidal", "pentagonal pyramidal", \
                    "pentagonal planar"]}
                                             
def get_molecular_formula():
    valid = False
    while not valid:
        formula = input("Molecular Formula (ex. CO2) : ")
        molecule = Molecule(formula, 0)
        molecule.get_atom_list(formula)
        if not molecule.atoms_are_nonmetals() or molecule.atom_list == []:
            print("Enter a molecular formula using only " +\
                  "nonmetals and their respective abundances " +\
                  "within the compound.")
        else:
            return formula

def get_charge():
    valid = False
    while not valid:
        charge = str(input("Molecule's Charge (ex. -1, 0, +2): "))
        if charge[0] == "+" or charge[0] == "-" or \
           charge[0] == "0":
            return charge

def is_finished():
    valid = False
    while not valid:
        finished = input("Are you finished? (yes/no): ")
        finished = finished.lower()
        if finished == "yes":
            return True
        elif finished == "no":
            return False
        else:
            print("Please enter either yes or no.")

def main():
    print("Welcome to the Lewis Diagram Program!", \
          "\nThis program is designed to accept a covalent compound", \
          "with one central atom, then display the molecule's Lewis Diagram.")
    finished = False
    while not finished:
        formula = get_molecular_formula()
        charge = get_charge()
        molecule = Molecule(formula, charge)
        molecule.get_atoms_bonds_lonepairs()
        # atom_bonds_lonepairs is a list containing 3 lists:
        # 1) the atoms in the molecule
        # 2) the number of bonds each atom has
        # 3) the number of lone pairs each atom has
        # The indexes in lists 2 & 3 correspond to their
        # respective atoms in list 1
        lewis = Lewis_Diagram(molecule.atoms_bonds_lonepairs,
                              molecule.central_atom_index)
        lewis.get_e_geo()
        lewis.get_mol_geo()
        lewis.draw_diagram(0, 0)
        print("Electron Geometry: " + lewis.e_geo)
        print("Molecular Geometry: " + lewis.mol_geo)
        finished = is_finished()
        s = turtle.Screen()
        s.clearscreen()
    print("Thank you for using the Lewis Diagram Program!")    
        
class Molecule:
    def __init__(self, formula, charge):
        self.formula = formula
        self.charge = charge
        self.atom_list = []
        self.atoms_bonds_lonepairs = []
        self.central_atom = ""
        self.central_atom_index = 0
        self.lonepairs_list = []

    def get_atoms_bonds_lonepairs(self):
        self.get_atom_list(self.formula)
        self.atoms_bonds_lonepairs.append(self.atom_list)
        self.get_central_atom()
        bonds = Bonds(self.charge, self.atom_list, \
                      self.central_atom_index)
        bonds.get_pos_neg_charge()
        bonds.get_e_shared()
        if self.is_diatomic():
            if self.has_hydrogen():
                bond_num = 1
            else:
                bond_num = bonds.get_bond_order()
            if self.central_atom_index == 0:
                self.atoms_bonds_lonepairs.append([0, bond_num])
            else:
                self.atoms_bonds_lonepairs.append([bond_num, 0])
            if self.has_hydrogen:
                self.h_diatomic_lonepairs()
            else:
                self.diatomic_lonepairs()
        else:
            bonds.get_bond_num()
            bonds.get_bonds_wanted_list()
            bonds.start_bonds_list()
            bonds.add_more_bonds()
            self.atoms_bonds_lonepairs.append(bonds.bonds_list)
            self.get_lonepairs_list(bonds.bonds_list)
        self.atoms_bonds_lonepairs.append(self.lonepairs_list)

    def get_atom_list(self, formula):
        i = 1
        while i <= len(formula):
            if formula[-i].isdigit():
                if formula[-(i + 1)].islower():
                    for n in range(
                        int(formula[- i])):
                        atom = formula[-(i + 2): - i]
                        self.atom_list.append(atom)                               
                    i += 3
                else:
                    for n in range(int(formula[- i])):
                        atom = formula[-(i + 1)]
                        self.atom_list.append(atom)
                    i += 2
            elif formula[-i].islower():
                self.atom_list.append(formula[-(i + 1)] + formula[-i])
                i += 2
            else:
                self.atom_list.append(formula[-i])
                i += 1
        return self.atom_list

    def atoms_are_nonmetals(self):
        for atom in self.atom_list:
            if atom not in atom_period_group_eneg:
                return False
        return True 

    def get_central_atom(self):
        # the central atom is the atom with the least electronegativity
        eneg_list = []
        for atom in self.atom_list:
            eneg_list.append(atom_period_group_eneg[atom][2])
        least_eneg = min(eneg_list)
        self.central_atom_index = eneg_list.index(least_eneg)
        self.central_atom = self.atom_list[self.central_atom_index]

    def get_lonepairs_list(self, bonds_list):
        bonds = Bonds(self.charge, self.atom_list, \
                      self.central_atom_index)
        bonds.get_pos_neg_charge()
        bonds.get_e_shared()
        for i in range(len(self.atom_list)):
            if i != self.central_atom_index:
                atom = Atom(self.atom_list[i], bonds_list[i])
                self.lonepairs_list.append(atom.lonepairs)
            else:
                if bonds.follows_octet_rule():
                    atom = Atom(self.atom_list[i], sum(bonds_list))
                    self.lonepairs_list.append(atom.lonepairs)
                else:
                    extra_e = (bonds.total_e_has - \
                                (8 * (len(self.atom_list) - 1)))
                    lonepairs = (extra_e // 2) + (extra_e % 2)
                    self.lonepairs_list.append(lonepairs)
            if self.lonepairs_list[i] < 0:
                self.lonepairs_list[i] = 0

    def is_diatomic(self):
        return len(self.atom_list) == 2

    def has_hydrogen(self):
        return "H" in self.atom_list

    def h_diatomic_lonepairs(self):
        self.lonepairs_list = []
        for i in range(len(self.atom_list)):
            if self.atom_list[i] == "H" or self.atom_list[i] == "He":
                self.lonepairs_list.append(0)
            else:
                self.lonepairs_list.append(3)
                # H and He's orbitals only hold 2 e-
                # other elements' orbitals hold 8 e-

    def diatomic_lonepairs(self):
        bonds = Bonds(self.charge, self.atom_list, \
                      self.central_atom_index)
        bonds.get_pos_neg_charge()
        bonds.get_e_shared()
        bonds.get_bond_order()
        lonepairs_total = (bonds.total_e_has // 2) - bonds.bond_order
        self.lonepairs_list = [lonepairs_total // 2, lonepairs_total // 2]
        if lonepairs_total % 2 == 1:
            for i in range(len(self.atom_list)):
                if i != self.central_atom_index:
                    self.lonepairs_list[i] += 1
    
class Bonds(Molecule):
    def __init__(self, charge, atom_list, central_atom_index):
        super().__init__(self, charge)
        self.total_e_has = 0
        self.total_e_full = 0
        self.e_shared = 0
        self.bonds_num = 0
        self.bond_order = 0
        # Due to the possiblity of having muliples of the same element,
        # a dictionary with elements as keys cannot be used.
        # Instead, I have used used lists whose indexes corrospond
        # to the atom_list. 
        self.bonds_wanted_list = []
        # atoms with a negative charge has extra electrons
        # atoms with a positive charge have fewer than normal electrons
        self.pos_charge = 0
        self.neg_charge = 0
        self.bonds_list = []
        self.atom_list = atom_list
        self.central_atom_index = central_atom_index

    def get_pos_neg_charge(self):
        if self.charge == 0:
            pass
        elif self.charge[0] == "+":
            self.pos_charge += int(self.charge[1:])
        elif self.charge[0] == "-":
            self.neg_charge += int(self.charge[1:])
        
    def get_e_shared(self):
        for element in self.atom_list:
            atom = Atom(element, 0)
            self.total_e_full += atom.e_full
            self.total_e_has += atom.e_has
        self.total_e_has += self.neg_charge
        self.total_e_has -= self.pos_charge
        self.e_shared += self.total_e_full - self.total_e_has

    def follows_octet_rule(self):
        return (len(self.atom_list) - 1) <= (self.e_shared / 2)

    def get_bond_num(self):
        if self.follows_octet_rule():
            self.bonds_num = self.e_shared // 2
        else:
            print("This molecule breaks the Octet Rule")
            self.bonds_num = len(self.atom_list) - 1

    def get_bonds_wanted_list(self):
        for element in self.atom_list:
            if element != self.atom_list[self.central_atom_index]:
                atom = Atom(element, 0)
                self.bonds_wanted_list.append(atom.e_needed)
            else:
                self.bonds_wanted_list.append(0)

    def start_bonds_list(self):
        for i in range(len(self.atom_list)):
            if i != self.central_atom_index:
                self.bonds_list.append(1)
                self.bonds_num -= 1
                self.bonds_wanted_list[i] -= 1
            else:
                self.bonds_list.append(0)
                self.bonds_wanted_list[
                    self.central_atom_index] = 0
                
    def add_more_bonds(self):
        while self.bonds_num > 0:
            max_need_index = self.bonds_wanted_list.index(
                max(self.bonds_wanted_list))
            self.bonds_list[max_need_index] += 1
            self.bonds_wanted_list[max_need_index] -= 1
            self.bonds_num -= 1

    def get_avg_atomic_num(self):
        atomic_num_total = 0
        for atom in self.atom_list:
            orbital = atom_period_group_eneg[atom][0]
            if  orbital == 1:
                atomic_num = atom_period_group_eneg[atom][1]
            elif orbital == 2:
                atomic_num = atom_period_group_eneg[atom][1] + orbital
            else:
                # the atomic number of elements above above the
                # second orbital level is arbitrary
                # because they all make the avg atomic num > 7
                atomic_num == 20
            atomic_num_total += atomic_num
        avg_atomic_num = (atomic_num_total // 2) + (atomic_num_total % 2)
        return avg_atomic_num

    def add_to_diagram(self, diagram, e, orbital, e_spots):
        for spot in range(e_spots):
            if e == 0:
                return 0
            else:
                diagram[orbital] += 1
                e -= 1
        return e

    def MO_diagram(self, diagram, orbital_list, e):
        for orbital in orbital_list:
            if e > 0:
                if orbital[-2:] == "pi":
                    e_spots = 4
                else:
                    e_spots = 2
                e = self.add_to_diagram(diagram, e, orbital, e_spots)
        return diagram
        
    def get_bond_order(self):
        e = self.total_e_has
        avg_atomic_num = self.get_avg_atomic_num()
        # MO_diagram stands for Molecular Orbital diagram
        # the order the MOs are filled depends on the avg_atomic_num
        # finding the e num in high and low energy MOs is needed
        # to calcualte the bond order, the number of
        # bonds between the elements
        start_diagram = {"low_sigma_s": 0, "high_sigma_s": 0, \
                         "low_pi": 0, "high_pi" : 0, \
                         "low_sigma_p": 0, "high_sigma_p": 0}
        if avg_atomic_num <= 7:
            orbital_list = ["low_sigma_s", "high_sigma_s", "low_pi", \
                            "low_sigma_p", "high_pi", "high_sigma_p"]
        else:
            orbital_list = ["low_sigma_s", "high_sigma_s", \
                            "low_sigma_p", "low_pi", \
                            "high_pi", "high_sigma_p"]
        full_diagram = self.MO_diagram(start_diagram, orbital_list, e)
        high_energy_e_total = full_diagram["high_sigma_s"] + \
                              full_diagram["high_pi"] + \
                              full_diagram["high_sigma_p"]
        low_energy_e_total = full_diagram["low_sigma_s"] + \
                             full_diagram["low_pi"] +\
                             full_diagram["low_sigma_p"]
        self.bond_order = ((low_energy_e_total - high_energy_e_total) // 2) +\
                          ((low_energy_e_total - high_energy_e_total) % 2)

class Atom:
    def __init__(self, element, bonds):
        if element == "H" or element == "He":
            self.e_full = 2
        else:
            self.e_full = 8
        self.e_has = atom_period_group_eneg[element][1]
        self.e_needed = self.e_full - self.e_has
        self.bonds = bonds
        self.extra_e = self.e_full - (self.bonds * 2)
        self.lonepairs = (self.extra_e // 2) + (self.extra_e % 2)

class Lewis_Diagram:
    def __init__(self, atoms_bonds_lonepairs, central_atom_index):
        self.mol_geometry = ""
        self.atoms_bonds_lonepairs = atoms_bonds_lonepairs
        self.atoms_list = self.atoms_bonds_lonepairs[0]
        self.bonds_list = self.atoms_bonds_lonepairs[1]
        self.lonepairs_list = self.atoms_bonds_lonepairs[2]
        self.central_atom_index = central_atom_index
        self.central_atom = self.atoms_list[self.central_atom_index]
        self.central_lonepairs = self.lonepairs_list[
            self.central_atom_index]
        self.e_geo = ""
        self.mol_geo = ""
        self.x = 0
        self.y = 0
        self.centralx = self.x + 15

    def get_e_geo(self):
        total_e_groups = (len(self.bonds_list) - 1) + \
                         self.central_lonepairs
        if len(self.atoms_list) == 2 or total_e_groups == 2:
            self.e_geo = "linear"
        elif total_e_groups == 3:
            self.e_geo = "trigonal planar"
        elif total_e_groups == 4:
            self.e_geo = "tetrahedral"
        elif total_e_groups == 5:
            self.e_geo = "trigonal bipyramidal"
        elif total_e_groups == 6:
            self.e_geo = "octahedral"
        elif total_e_groups == 7:
            self.e_geo = "pentagonal bipyramidal"

    def get_mol_geo(self):
        self.mol_geo = e_mol_geometry[self.e_geo][
            self.central_lonepairs]

    def linear_2atoms(self):
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        if self.central_lonepairs == 3:
            c.dotright(self.centralx, self.y)
            c.dotup(self.centralx, self.y)
            c.dotdown(self.centralx, self.y)
        elif self.central_lonepairs == 2:
            c.dotupright(self.centralx, self.y)
            c.dotdownright(self.centralx, self.y)
        elif self.central_lonepairs == 1:
            c.dotright(self.centralx, self.y)

    def linear_3atoms(self):
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        if self.central_lonepairs == 2:
            c.dotup(self.centralx, self.y)
            c.dotdown(self.centralx, self.y)
        elif self.central_lonepairs == 1:
            c.dotup(self.centralx, self.y)

    def draw_linear(self):
        n = 0
        for i in range(len(self.bonds_list)):
            if i != self.central_atom_index:
                n += 1
                d = Draw(self.x, self.y, self.atoms_list[i],
                         self.lonepairs_list[i])
                if n == 1:
                    if self.bonds_list[i] == 3:
                        d.left3()
                    elif self.bonds_list[i] == 2:
                        d.left2()
                    else:
                        d.left1()
                if n == 2:
                    if self.bonds_list[i] == 3:
                        d.right3()
                    elif self.bonds_list[i] == 2:
                        d.right2()
                    else:
                        d.right1()
        if n == 1:
            self.linear_2atoms()
        elif n == 2:
            self.linear_3atoms()

    def draw_trig_planar(self):
        n = 0
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        for i in range(len(self.bonds_list)):
            if i != self.central_atom_index:
                n += 1
                d = Draw(self.x, self.y, self.atoms_list[i],
                         self.lonepairs_list[i])
                if n == 1:
                    if self.bonds_list[i] == 2:
                        d.downleft2()
                    else:
                        d.downleft1()
                if n == 2:
                    if self.bonds_list[i] == 2:
                        d.downright2()
                    else:
                        d.downright1()
                if n == 3:
                    if self.bonds_list[i] == 2:
                        d.up2()
                    else:
                        d.up1()
            if self.central_lonepairs == 1:
                c.dotup(self.centralx, self.y)

    def draw_tetrahedral(self):
        n = 0
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        for i in range(len(self.bonds_list)):
            if i != self.central_atom_index:
                n += 1
                d = Draw(self.x, self.y, self.atoms_list[i],
                         self.lonepairs_list[i])
                if n == 1:
                    d.downright1()
                elif n == 2:
                    d.downleft1()
                elif n == 3:
                    d.down1()
                elif n == 4:
                    d.up1()
        if self.central_lonepairs == 2:
            c.dotupright(self.centralx, self.y)
            c.dotupleft(self.centralx, self.y)
        elif self.central_lonepairs == 1:
            c.dotup(self.centralx, self.y)

    def draw_trig_bipyramidal(self):
        n = 0
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        for i in range(len(self.bonds_list)):
            if i != self.central_atom_index:
                n += 1
                d = Draw(self.x, self.y, self.atoms_list[i],
                         self.lonepairs_list[i])
                if n == 1:
                    d.right1()
                elif n == 2:
                    d.left1()
                if self.mol_geo == "t-shaped":
                    if n == 3:
                       d.down1()
                elif self.mol_geo == "seesaw":
                    if n == 3:
                        d.downleft1()
                    elif n == 4:
                        d.downright1()
                elif self.mol_geo == self.e_geo:
                    if n == 3:
                        d.upleft1()
                    if n == 4:
                        d.down1()
                    elif n == 5:
                        d.up1()
        if self.central_lonepairs == 3:
            c.dotdownleft(self.centralx, self.y)
            c.dotdownright(self.centralx, self.y)
            c.dotup(self.centralx, self.y)
        elif self.central_lonepairs == 2:
            c.dotupleft(self.centralx, self.y)
            c.dotupright(self.centralx, self.y)
        elif self.central_lonepairs == 1:
            c.dotup(self.centralx, self.y)
                
    def draw_octahedral(self):
        n = 0
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        for i in range(len(self.bonds_list)):
            if i != self.central_atom_index:
                d = Draw(self.x, self.y, self.atoms_list[i],
                         self.lonepairs_list[i])
                n += 1
                if n == 1:
                    d.upright1()
                elif n == 2:
                    d.upleft1()
                elif n == 3:
                    d.downright1()
                elif n == 4:
                    d.downleft1()
                elif n == 5:
                    d.up1()
                elif n == 6:
                    d.down1()
        if self.central_lonepairs == 2:
            c.dotup(self.centralx, self.y)
            c.dotdown(self.centralx, self.y)
        elif self.central_lonepairs == 1:
            c.dotdown(self.centralx, self.y)

    def draw_pent_bipyramidal(self):
        n = 0
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        for i in range(len(self.bonds_list)):
            if i != self.central_atom_index:
                n += 1
                d = Draw(self.x, self.y, self.atoms_list[i],
                         self.lonepairs_list[i])
                if n == 1:
                    d.upright1()
                elif n == 2:
                    d.upleft1()
                elif n == 3:
                    d.downright1()
                elif n == 4:
                    d.downleft1()
                elif n == 5:
                    d.left1()
                elif n == 6:
                    d.down1()
                elif n == 7:
                    d.up1()
        if self.central_lonepairs >= 1:
            c.dotup(self.centralx, self.y)
        if self.central_lonepairs == 2:
            c.dotdown(self.centralx, self.y)

    def draw_diagram(self, x, y):
        x = self.x
        y = self.y
        c = Draw(self.centralx, self.y, self.central_atom,
                 self.central_lonepairs)
        c.central()
        if self.e_geo == "linear":
            self.draw_linear()
        elif self.e_geo == "trigonal planar":
            self.draw_trig_planar()
        elif self.e_geo == "tetrahedral":
            self.draw_tetrahedral()
        elif self.e_geo == "trigonal bipyramidal":
            self.draw_trig_bipyramidal()
        elif self.e_geo == "octahedral":
            self.draw_octahedral()
        elif self.e_geo == "pentagonal bipyramidal":
            self.draw_pent_bipyramidal()

class Draw:
    def __init__(self, x, y, atom, lonepairs):
        self.x = x
        self.xa = x + 22
        self.xb = x + 25
        self.xc = x + 28
        self.y = y
        self.ya = x + 12
        self.yb = x + 15
        self.yc = x + 18
        self.atom = atom
        self.atom_len = len(atom)
        self.lonepairs = lonepairs
        self.t = turtle.Turtle()

    def central(self):
        self.t.penup()
        self.t.goto(self.x, self.y)
        self.t.write(self.atom, font = ("Arial", 20, "normal"))
        self.t.ht()

    def draw1(self, x1, y1, x2, y2, xwrite, ywrite):
        self.t.penup()
        self.t.goto(x1, y1)
        self.t.pendown()
        self.t.goto(x2, y2)
        self.t.penup()
        self.t.goto(xwrite, ywrite)
        self.t.pendown()
        self.t.write(self.atom, font = ("Arial", 20, "normal"))
        self.t.ht()

    def right1(self):
        x1 = self.x + 50
        x2 = self.x + 80
        self.draw1(x1, self.yb, x2, self.yb, x2, self.y)
        if self.lonepairs == 1:
            self.dotright(x2, self.y)
        elif self.lonepairs == 2:
            self.dotupright(x2, self.y)
            self.dotdownright(x2, self.y)
        elif self.lonepairs == 3:
            self.dotup(x2, self.y)
            self.dotright(x2, self.y)
            self.dotdown(x2, self.y)

    def left1(self):
        x1 = self.x - 5
        x2 = self.x - 35
        if self.atom_len == 1:
            xwrite = x2 - 20
        elif self.atom_len == 2:
            xwrite = x2 - 40
        self.draw1(x1, self.yb, x2, self.yb, xwrite, self.y)
        if self.lonepairs == 1:
            self.dotleft(xwrite, self.y)
        elif self.lonepairs == 2:
            self.dotupleft(xwrite, self.y)
            self.dotdownleft(xwrite, self.y)
        elif self.lonepairs == 3:
            self.dotup(xwrite, self.y)
            self.dotleft(xwrite, self.y)
            self.dotdown(xwrite, self.y)

    def up1(self):
        y1 = self.y + 35
        y2 = self.y + 65
        self.draw1(self.xb, y1, self.xb, y2, self.xa, y2)
        if self.lonepairs == 1:
            self.dotup(self.xa, y2)
        elif self.lonepairs == 2:
            self.dotupleft(self.xa, y2)
            self.dotupright(self.xa, y2)
        elif self.lonepairs == 3:
            self.dotup(self.xa, y2)
            self.dotleft(self.xa, y2)
            self.dotright(self.xa, y2)

    def down1(self):
        y2 = self.y - 30
        self.draw1(self.xb, self.y, self.xb, y2, self.xa, y2 - 25)
        if self.lonepairs == 1:
            self.dotdown(self.xa, y2 - 25)
        elif self.lonepairs == 2:
            self.dotupleft(self.xa, y2 - 25)
            self.dotupright(self.xa, y2 - 25)
        elif self.lonepairs == 3:
            self.dotdown(self.xa, y2 - 25)
            self.dotleft(self.xa, y2 - 25)
            self.dotright(self.xa, y2 - 25)

    def upright1(self):
        x1 = self.x + 50
        x2 = self.x + 80
        y1 = self.y + 35
        y2 = self.y + 65
        self.draw1(x1, y1, x2, y2, x2, y2) 
        if self.lonepairs == 1:
            self.dotright(x2, y2)
        elif self.lonepairs == 2:
            self.dotup(x2, y2)
            self.dotright(x2, y2)
        elif self.lonepairs == 3:
            self.dotup(x2, y2)
            self.dotright(x2, y2)
            self.dotleft(x2, y2)
            
    def upleft1(self):
        x1 = self.x - 5
        x2 = self.x - 35
        y1 = self.y + 35
        y2 = self.y + 65
        if self.atom_len == 1:
            xwrite = x2 - 20
        elif self.atom_len == 2:
            xwrite = x2 - 40
        self.draw1(x1, y1, x2, y2, xwrite, y2)
        if self.lonepairs == 1:
            self.dotleft(xwrite, y2)
        elif self.lonepairs == 2:
            self.dotup(xwrite, y2)
            self.dotleft(xwrite, y2)
        elif self.lonepairs == 3:
            self.dotup(xwrite, y2)
            self.dotright(xwrite, y2)
            self.dotleft(xwrite, y2)

    def downright1(self):
        x1 = self.x + 50
        x2 = self.x + 80
        y2 = self.y - 30
        self.draw1(x1, self.y, x2, y2, x2, y2 - 25)
        if self.lonepairs == 1:
            self.dotright(x2, y2 - 25)
        elif self.lonepairs == 2:
            self.dotdown(x2, y2 - 25)
            self.dotright(self.xa, y2 - 25)
        elif self.lonepairs == 3:
            self.dotdown(x2, y2 - 25)
            self.dotleft(x2, y2 - 25)
            self.dotright(x2, y2 - 25)

    def downleft1(self):
        x1 = self.x - 5
        x2 = self.x - 35
        y2 = self.y - 30
        if self.atom_len == 1:
            xwrite = x2 - 20
        elif self.atom_len == 2:
            xwrite = x2 - 40
        self.draw1(x1, self.y, x2, y2, xwrite, y2 - 25)
        if self.lonepairs == 1:
            self.dotleft(xwrite, y2 - 25)
        elif self.lonepairs == 2:
            self.dotdown(xwrite, y2 - 25)
            self.dotleft(xwrite, y2 - 25)
        elif self.lonepairs == 3:
            self.dotdown(xwrite, y2 - 25)
            self.dotleft(xwrite, y2 - 25)
            self.dotright(xwrite, y2 - 25)

    def draw2(self, x1a, x1b, y1a, y1b, x2a, x2b,
              y2a, y2b, xwrite, ywrite):
        self.t.penup()
        self.t.goto(x1a, y1a)
        self.t.pendown()
        self.t.goto(x2a, y2a)
        self.t.penup()
        self.t.goto(x1b, y1b)
        self.t.pendown()
        self.t.goto(x2b, y2b)
        self.t.penup()
        self.t.goto(xwrite, ywrite)
        self.t.pendown()
        self.t.write(self.atom, font = ("Arial", 20, "normal"))
        self.t.ht()

    def right2(self):
        x1 = self.x + 50
        x2 = self.x + 80
        self.draw2(x1, x1, self.ya, self.yc, x2, x2,
                   self.ya, self.yc, x2, self.y)
        if self.lonepairs == 1:
            self.dotright(x2, self.y)
        elif self.lonepairs == 2:
            self.dotupright(x2, self.y)
            self.dotdownright(x2, self.y)

    def left2(self):
        x1 = self.x - 5
        x2 = self.x - 35
        if self.atom_len == 1:
            xwrite = x2 - 20
        elif self.atom_len == 2:
            xwrite = x2 - 40
        self.draw2(x1, x1, self.ya, self.yc, x2, x2,
                   self.ya, self.yc, xwrite, self.y)
        if self.lonepairs == 1:
            self.dotleft(xwrite, self.y)
        elif self.lonepairs == 2:
            self.dotupleft(xwrite, self.y)
            self.dotdownleft(xwrite, self.y)
        elif self.lonepairs == 3:
            self.dotdown(xwrite, self.y)
            self.dotup(xwrite, self.y)
            self.dotleft(xwrite, self.y)

    def up2(self):
        y1 = self.y + 35
        y2 = self.y + 65
        self.draw2(self.xa, self.xc, y1, y1, self.xa,
                   self.xc, y2, y2, self.xa, y2)
        if self.lonepairs == 1:
            self.dotup(self.xa, y2)
        elif self.lonepairs == 2:
            self.dotupleft(self.xa, y2)
            self.dotupright(self.xa, y2)

    def down2(self):
        y1 = self.y - 35
        y2 = self.y - 65
        self.draw2(self.xa, self.xc, y1, y1, self.xa,
                   self.xc, y2, y2, self.xa, y2) 
        if self.lonepairs == 1:
            self.dotdown(self.xa, y2 - 25)
        elif self.lonepairs == 2:
            self.dotupleft(self.xa, y2 - 25)
            self.dotupright(self.xa, y2 - 25)

    def downright2(self):
        x1a = self.x + 47
        x1b = self.x + 53
        x2a = self.x + 77
        x2b = self.x + 83
        y2 = self.y - 30
        self.draw2(x1a, x1b, self.y, self.y, x2a, x2b, y2, \
                  y2, x2a, y2 - 25)
        if self.lonepairs == 1:
            self.dotright(x2a, y2 - 25)
        elif self.lonepairs == 2:
            self.dotdown(x2a, y2 - 25)
            self.dotright(x2a, y2 - 25)

    def downleft2(self):
        x1a = self.x - 2
        x1b = self.x - 8
        x2a = self.x - 32
        x2b = self.x - 38
        y2 = self.y - 30
        if self.atom_len == 1:
            xwrite = x2a - 20
        elif self.atom_len == 2:
            xwrite = x2a - 40
        self.draw2(x1a, x1b, self.y, self.y, x2a, x2b, \
                   y2, y2, xwrite, y2 - 25)
        if self.lonepairs == 1:
            self.dotleft(xwrite, y2 - 25)
        elif self.lonepairs == 2:
            self.dotdown(xwrite, y2 - 25)
            self.dotleft(xwrite, y2 - 25)

    def draw3(self, x1a, x1b, x1c, y1a, y1b, \
              y1c, x2a, x2b, x2c, y2a, y2b, y2c, \
              xwrite, ywrite):
        self.draw2(x1a, x1b, y1a, y1b, x2a, x2b, \
                   y2a, y2b, xwrite, ywrite)
        self.draw1(x1c, y1c, x2c, y2c, xwrite, ywrite)
        self.t.penup()
        self.t.goto(xwrite, ywrite)
        self.t.pendown()
        self.t.write(self.atom, font = ("Arial", 20, "normal"))

    def right3(self):
        x1 = self.x + 50
        x2 = self.x + 80
        self.draw3(x1, x1, x1, self.ya, self.yb, \
                   self.yc, x2, x2, x2, self.ya, self.yb, \
                   self.yc, x2, self.y)
        if self.lonepairs == 1:
            self.dotright(x2, self.y)

    def left3(self):
        x1 = self.x - 5
        x2 = self.x - 35
        if self.atom_len == 1:
            xwrite = x2 - 20
        elif self.atom_len == 2:
            xwrite = x2 - 40
        self.draw3(x1, x1, x1, self.ya, self.yb, \
                   self.yc, x2, x2, x2, self.ya, self.yb, \
                   self.yc, xwrite, self.y)
        if self.lonepairs == 1:
            self.dotleft(xwrite, self.y)

    def dot(self, x1, y1, x2, y2):
        self.t.penup()
        self.t.goto(x1, y1)
        self.t.dot()
        self.t.goto(x2, y2)
        self.t.dot()
        self.t.ht()

    def dotup(self, x, y):
        self.dot(x + 2, y + 33, x + 12, y + 33)

    def dotdown(self, x, y):
        self.dot(x + 2, y, x + 12, y)

    def dotleft(self, x, y):
        self.dot(x - 5, y + 10, x - 5, y + 20)

    def dotright(self, x, y):
        self.dot(x + 27, y + 10, x + 27, y + 20)

    def dotupright(self, x, y):
        self.dot(x + 17, y + 32, x + 23, y + 30)

    def dotupleft(self, x, y):
        self.dot(x - 3, y + 30, x + 2, y + 32)

    def dotdownright(self, x, y):
        self.dot(x + 17, y, x + 23, y + 5)

    def dotdownleft(self, x, y):
        self.dot(x - 4, y + 5, x + 2, y)

main()
