"""
This is a set of useful features for manipulating basic lammps output files

###############################################################################
###############################################################################
author: Aaron Schwan
email: schwanaaron@gmail.com
github: https://github.com/AaronSchwan
website: https://aaronschwan.wordpress.com/lammps-file-manipulation-package/
###############################################################################
###############################################################################

"""


#default imports
import sys
import os
import time
import warnings
import types
import copy

#non-default imports
import pandas as pd
import numpy as np

################################################################################
#Dealing with lammps dump files#################################################
################################################################################

class dumpFile:

    """
    This is a class optimized for operations on dump files exported from lammps
    programs

    obj = dumpFile(timestep:int,numberofatoms:int,boxbounds:pd.DataFrame,atoms:pd.DataFrame,serial=None)

    Alternative class construction methods(file_path = path to file):
        dumpFile.lammps_dump(cls, file_path) #reads in a standard lammps dump
        dumpFile.pandas_to_dumpfile(cls, file_path) #reads in lammps data from pandas dataframe to new class

    valid property calls:
        obj.timestep = returns timestep in the file[int]
        obj.numberofatoms = numbers of atoms in the dump[int]
        obj.atoms = atomic data[pd.DataFrame]
        obj.boxbounds = returns the bounds with type,low,high in a pandas datframe[pd.DataFrame]
        obj.volume = returns cubic shaped volume
        obj.xy_area = returns square shaped area of xy plane
        obj.xz_area = returns square shaped area of xz plane
        obj.yz_area = returns square shaped area of yz plane
        obj.xlo = low x value in obj.atoms
        obj.ylo = low y value in obj.atoms
        obj.zlo = low z value in obj.atoms
        obj.xhi = high x value in obj.atoms
        obj.yhi = high y value in obj.atoms
        obj.zhi = high z value in obj.atoms



    method calls:
        obj.translate() = method to translate atom locations of the instance this
                            will return a new instance of the class in order to
                            preserve data integrety for the class instance

    dunder calls:
        "obj1 == obj2" = returns if the atomic positional distances are identical
                            uses the class variable checking_tolerance for amount
                            of precission in check

        "merge_obj = obj1+obj2" = alternative way to merge data of the obj.atoms
    """

    #precision based variables
    class_tolerance = 12 #the accuarcy of the classes operational functions
    checking_tolerance = 3 #how many decimals the classes attributes will be checked to

    ##atomic identification
    id = "id"
    type = "type"
    ##cartesian
    x_axis_cart = "x"
    y_axis_cart = "y"
    z_axis_cart = "z"


    def __init__(self,timestep:int,boundingtypes:list,atoms:pd.DataFrame):
        self.timestep = timestep
        self.boundingtypes = boundingtypes
        self.atoms = atoms

    #property defined functions#################################################
    @property
    def numberofatoms(self):
        #get axis
        id = self.id
        #finds the number of atoms for the given data
        return len(self.atoms[id])

    @property
    def boxbounds(self):
        #get axes
        x = self.x_axis_cart
        y = self.y_axis_cart
        z = self.z_axis_cart

        #getting bounds and making custom format to values
        box_dict = {"labels":["type","low","high"],"x":[self.boundingtypes[0],min(self.atoms[x]),max(self.atoms[x])],"y":[self.boundingtypes[1],min(self.atoms[y]),max(self.atoms[y])],"z":[self.boundingtypes[2],min(self.atoms[z]),max(self.atoms[z])]}
        return pd.DataFrame.from_dict(box_dict).set_index("labels")

    @property
    def xlo(self):
        #get axis
        x = self.x_axis_cart
        return min(self.atoms[x])

    @property
    def xhi(self):
        #get axis
        x = self.x_axis_cart
        return max(self.atoms[x])

    @property
    def ylo(self):
        #get axis
        y = self.y_axis_cart
        return min(self.atoms[y])

    @property
    def yhi(self):
        #get axis
        y = self.y_axis_cart
        return max(self.atoms[y])

    @property
    def zlo(self):
        #get axis
        z = self.z_axis_cart
        return min(self.atoms[z])

    @property
    def zhi(self):
        #get axis
        z = self.z_axis_cart
        return max(self.atoms[z])

    @property
    def volume(self):
        """
        gets the volume of the overall simulation cell as a box
        """
        x_range = abs(self.boxbounds.loc["high","x"] - self.boxbounds.loc["low","x"])
        y_range = abs(self.boxbounds.loc["high","y"] - self.boxbounds.loc["low","y"])
        z_range = abs(self.boxbounds.loc["high","z"] - self.boxbounds.loc["low","z"])

        return x_range*y_range*z_range

    @property
    def xy_area(self):
        """
        gets the area of the xy plane of the simulation cell as a box
        """
        x_range = abs(self.boxbounds.loc["high","x"] - self.boxbounds.loc["low","x"])
        y_range = abs(self.boxbounds.loc["high","y"] - self.boxbounds.loc["low","y"])

        return x_range*y_range

    @property
    def xz_area(self):
        """
        gets the area of the xy plane of the simulation cell as a box
        """
        x_range = abs(self.boxbounds.loc["high","x"] - self.boxbounds.loc["low","x"])
        z_range = abs(self.boxbounds.loc["high","z"] - self.boxbounds.loc["low","z"])

        return x_range*z_range

    @property
    def yz_area(self):
        """
        gets the area of the xy plane of the simulation cell as a box
        """
        z_range = abs(self.boxbounds.loc["high","z"] - self.boxbounds.loc["low","z"])
        y_range = abs(self.boxbounds.loc["high","y"] - self.boxbounds.loc["low","y"])

        return z_range*y_range



    #Alternative CLass Constructive Methods#####################################
    @classmethod
    def lammps_dump(cls,file_path:str):
        """
        uses path of raw lammps file **Must be a singular timestep

        will create the class of dumpFile once processed
        """
        raw_data = pd.read_csv(file_path,header = None)#getting data
        indexes = raw_data.index[raw_data[0].str.contains("ITEM: TIMESTEP")].tolist()#allowing check for singular

        if len(indexes) == 1:
            titles = raw_data.iloc[8].str.split(expand = True).iloc[0][2:].tolist()#getting titles of atomic data
            data = pd.DataFrame(raw_data.iloc[9:,0].str.split(' ',len(titles)-1).tolist(), columns = titles)#getting atomic data
            data =  data.apply(pd.to_numeric)

            #getting bounds and making custom format to values
            boxboundtype =  raw_data.iloc[4,0].replace("ITEM: BOX BOUNDS ","").split(" ")

            #returning class
            return cls(int(raw_data.iloc[1,0]),boxboundtype,data)

        elif len(indexes) > 1:
            raise Exception("FILE IMPORT ERROR: You may not import a multiple timestep file using this method please use the multiple_timestep_singular_file_dumps function")

        else:
            raise Exception("FILE IMPORT ERROR: check file formatting ")

    @classmethod
    def pandas_to_dumpfile(cls,raw_data:pd.DataFrame):
        """
        takes in a lammps dump file in the form of a singular column singular time step

        **Must include all data from first row "ITEM: TIMESTEP" to last row in one column
        """

        indexes = raw_data.index[raw_data[0].str.contains("ITEM: TIMESTEP")].tolist()#allowing check for singular

        if len(indexes) == 1:
            titles = raw_data.iloc[8].str.split(expand = True).iloc[0][2:].tolist()#getting titles of atomic data
            data = pd.DataFrame(raw_data.iloc[9:,0].str.split(' ',len(titles)-1).tolist(), columns = titles)#getting atomic data
            data =  data.apply(pd.to_numeric)

            #getting bounds and making custom format to values
            boxboundtype =  raw_data.iloc[4,0].replace("ITEM: BOX BOUNDS ","").split(" ")

            #returning class
            return cls(int(raw_data.iloc[1,0]),boxboundtype,data)

        elif len(indexes) > 1:
            raise Exception("FILE IMPORT ERROR: You may not import a multiple timestep file using this method please use the multiple_timestep_singular_file_dumps function")

        else:
            raise Exception("FILE IMPORT ERROR: check file formatting ")

    #Class methods##############################################################
    @classmethod
    def change_checking_tolerance(cls,value):
        cls.checking_tolerance = value

    @classmethod
    def change_class_tolerance(cls,value):
        cls.class_tolerance = value


    #dubble under functions#####################################################
    def __repr__(self):
        #returning the atoms by default when calling the function alone
        return "{TimeStep:"+str(self.timestep)+"\nBoundings"+str(self.boxbounds)+"\nColumns of atomic data"+str(self.atoms.columns)+"}"

    def __eq__(self,other):
        """
        This checks if the base conditions are equal such as the number of atoms,
        the boundary after transform to the first boundary, then if it passes
        both of those conditions it will check the atoms positions within a given
        tolerance(class variable name = class_tolerance)

        This will auto transform both dumpFile class instances to the positive
        quadrent system in order to avoid issues involving one transformed and
        one non transformed

        The comparison will be done using the leftmost class instances
        active coordinate system

        checks in order: number of atoms -> atomic positions
        """

        if self.numberofatoms == other.numberofatoms:
            #subtraction method because pd.equals will return false unless more operations are done
            df_self = self.translate(1).atoms[[self.id,self.x_axis_cart,self.y_axis_cart,self.z_axis_cart]]
            df_other = other.translate(1).atoms[[self.id,self.x_axis_cart,self.y_axis_cart,self.z_axis_cart]]
            df = (df_self-df_other).round(dumpFile.checking_tolerance)
            if df[self.id].eq(0).all() and df[self.x_axis_cart].eq(0).all()and df[self.y_axis_cart].eq(0).all() and df[self.z_axis_cart].eq(0).all():
                return True
            else:
                return False

        else:
            return False


    def __add__(self,other):
        """
        This is an alternative merge method first the class checks the compatability
        of the merge

        This acts be leaving the columns in the leftmost class instance untouched
        and appends the unique data to the leftmost class instance returning a new
        instance
        """
        if self == other:
            unique_columns = np.setdiff1d(other.atoms.columns.tolist(),self.atoms.columns.tolist())

            atomic_data = self.atoms.join(other.atoms[unique_columns])#merged atoms
            return dumpFile(self.timestep,self.boundingtypes,atomic_data)

        else:
            raise Exception("You may not add two classes where the atomic conditions/placements are not equal")

    #Class functional methods###################################################
    def translate(self,translation_operation,coordinate_sys = "active"):
        """
        This function transforms the atoms of the class to different quadrents
        labeled below.

        proper call:
        translated = class_instance.translate(translation_operation)

        The predefined quadrent operations will make sure one boundry is a 0 0 0

        Standards Cartesian:
         quadrent = x y z
                0 = centered at 0 0 0
                1 = + + +
                2 = + + -
                3 = + - +
                4 = + - -
                5 = - + +
                6 = - + -
                7 = - - +
                8 = - - -

        Custum Transform:
         The value is added to the atoms direction from the list in order [x_shift, y_shift, z_shift]

        """


        #defining new object
        dump_class_object = copy.deepcopy(self)


        id = self.id
        x = self.x_axis_cart
        y = self.y_axis_cart
        z = self.z_axis_cart
        #getting correction direction for minimums
        if min(dump_class_object.atoms[x]) > 0:
            dir_x = -1
        else:
            dir_x = 1

        if min(dump_class_object.atoms[y]) > 0:
            dir_y = -1
        else:
            dir_y = 1

        if min(dump_class_object.atoms[z]) > 0:
            dir_z = -1
        else:
            dir_z = 1

        #standard translation_operationtransforms
        if translation_operation == 1:

            dump_class_object.atoms[x] = dump_class_object.atoms[x] - min(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - min(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - min(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 2:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - min(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - min(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] + dir_z*min(dump_class_object.atoms[z])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - max(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 3:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - min(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] + dir_x*min(dump_class_object.atoms[y])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - max(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - min(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 4:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - min(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] + dir_x*min(dump_class_object.atoms[y])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - max(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] + dir_x*min(dump_class_object.atoms[z])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - max(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 5:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] + dir_x*min(dump_class_object.atoms[x])
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - max(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - min(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - min(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 6:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] + dir_x*min(dump_class_object.atoms[x])
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - max(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - min(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] + dir_x*min(dump_class_object.atoms[z])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - max(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 7:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] + dir_x*min(dump_class_object.atoms[x])
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - max(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] + dir_x*min(dump_class_object.atoms[y])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - max(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - min(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 8:
            dump_class_object.atoms[x] = dump_class_object.atoms[x] + dir_x*min(dump_class_object.atoms[x])
            dump_class_object.atoms[x] = dump_class_object.atoms[x] - max(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] + dir_x*min(dump_class_object.atoms[y])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - max(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] + dir_x*min(dump_class_object.atoms[z])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - max(dump_class_object.atoms[z])

            return dump_class_object

        elif translation_operation == 0:

            dump_class_object.atoms[x] = dump_class_object.atoms[x] - min(dump_class_object.atoms[x])
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - min(dump_class_object.atoms[y])
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - min(dump_class_object.atoms[z])

            dump_class_object.atoms[x] = dump_class_object.atoms[x] - (max(dump_class_object.atoms[x])/2)
            dump_class_object.atoms[y] = dump_class_object.atoms[y] - (max(dump_class_object.atoms[y])/2)
            dump_class_object.atoms[z] = dump_class_object.atoms[z] - (max(dump_class_object.atoms[z])/2)

            return dump_class_object

        elif type(translation_operation) == list:
            if all(isinstance(i, (float, int)) for i in translation_operation) and len(translation_operation) == 3:
                dump_class_object.atoms[x] = dump_class_object.atoms[x] + translation_operation[0]
                dump_class_object.atoms[y] = dump_class_object.atoms[y] + translation_operation[1]
                dump_class_object.atoms[z] = dump_class_object.atoms[z] + translation_operation[2]

                return dump_class_object


            else:
                 raise Exception("Not a valid input to translation function custom list")

        else:
             raise Exception("Not a valid input to translation function")

    #writing out functions
    def write_dump_file(self,file_path:str,mode:str = "a"):
        """
        This takes in a file path and writes a dumpFile class to the file path in
        standard lammps format

        write_lammps_dump(file_path:str,dump_class:dumpFile,mode:str = "a")

        file_path = path to file [str]
        mode = overwrite("w") or append("a") **default append [str]

        """

        if mode == "a" or mode == "w":
            precision = dumpFile.class_tolerance#get writing precision
            #defining new object
            dump_class_object = copy.deepcopy(self)

            with open(file_path,mode) as file:
                file.write("ITEM: TIMESTEP \n")
                file.write(str(dump_class_object.timestep))
                file.write("\n")
                file.write("ITEM: NUMBER OF ATOMS \n")
                file.write(str(dump_class_object.numberofatoms))
                file.write("\n")
                file.write("ITEM: BOX BOUNDS ")
                file.write(str(dump_class_object.boundingtypes[0])+" "+str(dump_class_object.boundingtypes[1])+" "+str(dump_class_object.boundingtypes[2])+"\n")#bound types
                file.write(str(round(dump_class_object.xlo,precision))+" "+str(round(dump_class_object.xhi,precision))+"\n")
                file.write(str(round(dump_class_object.ylo,precision))+" "+str(round(dump_class_object.yhi,precision))+"\n")
                file.write(str(round(dump_class_object.zlo,precision))+" "+str(round(dump_class_object.zhi,precision))+"\n")
                file.write("ITEM: ATOMS ")


            dump_class_object.atoms.round(precision).to_csv(file_path,mode = "a", index = False,sep = ' ')
            del dump_class_object


        else:
            raise Exception('Mode entered for writing is not recognized ["a"= append to files, "w"= overwrite file]')


    def write_dump_to_data_format(self, file_path:str,mode:str = "a"):
        """
        writes dumpFile class to a data file format

        mode = overwrite("w") or append("a") **default append [str]

        **this will only save the positions in data fromatting
        **primary use to write a initiallization data file for a lammps

        ** example

        # LAMMPS data file written by LammpsFileManipulation.py
        275184 atoms
        2 atom types
        0.4892064609 119.789657019 xlo xhi
        -158.7268972078 158.7268972078 ylo yhi
        0.4917072972 119.7871561826 zlo zhi

        Atoms  # atomic

        1 2 2.59911 -158.671 2.89486
        .
        .
        .

        """
        if mode == "a" or mode == "w":
            precision = dumpFile.class_tolerance#get writing precision
            #defining new object
            dump_class_object = copy.deepcopy(self)

            id = dump_class_object.id
            type = dump_class_object.type
            x = dump_class_object.x_axis_cart
            y = dump_class_object.y_axis_cart
            z = dump_class_object.z_axis_cart

            with open(file_path,"w") as file:
                file.write("# LAMMPS data file written by LammpsFileManipulation.py \n")
                file.write(str(dump_class.numberofatoms))
                file.write(" atoms \n")
                file.write(str(max(dump_class.atoms["type"])))
                file.write(" atom types \n")
                file.write(str(round(dump_class_object.xlo,precision))+" "+str(round(dump_class_object.xhi,precision))+" xlo xhi\n")
                file.write(str(round(dump_class_object.ylo,precision))+" "+str(round(dump_class_object.yhi,precision))+" ylo yhi\n")
                file.write(str(round(dump_class_object.zlo,precision))+" "+str(round(dump_class_object.zhi,precision))+" zlo zhi\n")
                file.write("\n\n")
                file.write("Atoms  # atomic\n\n")

            dump_class.atoms[[id,type, x, y, z]].round(precision).to_csv(file_path,mode = "a", index = False,header = False ,sep = ' ')
            del dump_class_object

        else:
            raise Exception('Mode entered for writing is not recognized ["a"= append to files, "w"= overwrite file]')

def group_translate(dump_files, translation_operation):
    """
    This takes in a group of dumps in the dictionary format of class and ####translates
    them as a group the same amount###

    because of the nature of this operation the values must be found first therefore
    when not using a custom transform the function loops through adding a small
    amount of time to this practice

    proper call:
    translated = group_translate(dump_files,quadrent)

    dump_files = {id:class,...}

    The predefined quadrent operations will make sure one boundry is a 0 0 0

    quadrent = x y z
        0 = centered at 0 0 0
        1 = + + +
        2 = + + -
        3 = + - +
        4 = + - -
        5 = - + +
        6 = - + -
        7 = - - +
        8 = - - -

    Custum Transform:
     The value is added to the atoms direction from the list in order [x_shift, y_shift, z_shift]

    """
    #if custom translation

    #getting max and min boundings
    id = dumpFile.id
    x = dumpFile.x_axis_cart
    y = dumpFile.y_axis_cart
    z = dumpFile.z_axis_cart

    if type(translation_operation) == list:
        if all(isinstance(i, (float, int)) for i in translation_operation) and len(translation_operation) == 3:
            tran_x = translation_operation[0]
            tran_y = translation_operation[1]
            tran_z = translation_operation[2]

        else:
             raise Exception("Not a valid input to translation function custom list")
    else:

        for ind,dump_class_id in enumerate(dump_files):
            dump_class = dump_files[dump_class_id]
            atomic_data = dump_class.atoms
            #get comparision values
            x_low_t = min(atomic_data[x])
            y_low_t = min(atomic_data[y])
            z_low_t = min(atomic_data[z])
            x_max_t = max(atomic_data[x])
            y_max_t = max(atomic_data[y])
            z_max_t = max(atomic_data[z])
            x_sh_max_t = max(atomic_data[x])-min(atomic_data[x])
            y_sh_max_t = max(atomic_data[y])-min(atomic_data[y])
            z_sh_max_t = max(atomic_data[z])-min(atomic_data[z])

            if ind == 0:
                #establishes starting point
                x_low = x_low_t
                y_low = y_low_t
                z_low = z_low_t
                x_max = x_max_t
                y_max = y_max_t
                z_max = z_max_t
                x_sh_max = x_sh_max_t
                y_sh_max = y_sh_max_t
                z_sh_max = z_sh_max_t


            #updating values if needed
            if x_low_t < x_low:
                x_low = x_low_t
            if y_low_t < y_low:
                y_low = y_low_t
            if z_low_t < z_low:
                z_low = z_low_t
            if x_max_t > x_max:
                x_max = x_max_t
            if y_max_t > y_max:
                y_max = y_max_t
            if z_max_t > z_max:
                z_max = z_max_t
            if x_sh_max_t > x_sh_max:
                x_sh_max = x_sh_max_t
            if y_sh_max_t > y_sh_max:
                y_sh_max = y_sh_max_t
            if z_sh_max_t > z_sh_max:
                z_sh_max = z_sh_max_t

        #getting correction direction for minimums
        if x_low > 0:
            dir_x = -1
        else:
            dir_x = 1

        if y_low > 0:
            dir_y = -1
        else:
            dir_y = 1

        if z_low > 0:
            dir_z = -1
        else:
            dir_z = 1

        #translating

        if translation_operation == 0:
            tran_x = -(x_max-x_low)/2 + dir_x*x_low
            tran_y = -(y_max-y_low)/2 + dir_y*y_low
            tran_z = -(z_max-z_low)/2 + dir_z*z_low

        elif translation_operation == 1:
            tran_x = -x_low
            tran_y = -y_low
            tran_z = -z_low

        elif translation_operation == 2:
            tran_x = -x_low
            tran_y = -y_low
            tran_z = -z_low-z_sh_max

        elif translation_operation == 3:
            tran_x = -x_low
            tran_y = -y_low-y_sh_max
            tran_z = -z_low

        elif translation_operation == 4:
            tran_x = -x_low
            tran_y = -y_low-y_sh_max
            tran_z = -z_low-z_sh_max

        elif translation_operation == 5:
            tran_x =  -x_low-x_sh_max
            tran_y = -y_low
            tran_z = -z_low

        elif translation_operation == 6:
            tran_x =  -x_low-x_sh_max
            tran_y = -y_low
            tran_z = -z_low-z_sh_max

        elif translation_operation == 7:
            tran_x =  -x_low-x_sh_max
            tran_y = -y_low-y_sh_max
            tran_z = -z_low

        elif translation_operation == 8:
            tran_x = -x_low-x_sh_max
            tran_y = -y_low-y_sh_max
            tran_z = -z_low-z_sh_max



        else:
             raise Exception("Not a valid input to translation function")

    #transforming classes
    translated_dump_files = {}
    for dump_class_id in dump_files:
        translated_dump_files[dump_class_id] = dump_files[dump_class_id].translate([tran_x,tran_y,tran_z])

    return translated_dump_files


def multiple_timestep_singular_file_dumps(file_path:str,ids:list = ["TimestepDefault"]):
    """
    this opens a multi-timestep lammps dump and converts it to a dictionary of
    dumpFile classes with the keys set to the timesteps

    ids:list = ["TimestepDefault"]
    ids are set to the dumpclass timestep by default however if there are duplicates
    this will override the timesteps so you can define the ids for the dictionary
    """
    dump_files = {} #dictionary of class

    raw_data = pd.read_csv(file_path,header = None)#getting data
    indexes = raw_data.index[raw_data[0] == "ITEM: TIMESTEP"].tolist()#getting splitting indexes

    if len(ids) == len(indexes) or ids == ["TimestepDefault"]:

        #splitting and iterating through pandas dataFrame
        for ind, index in enumerate(indexes):

            if ind < len(indexes)-1:
                #all except last index
                df = raw_data.loc[index:indexes[ind+1]-1,:]#making new dataFrame
                dump_class = dumpFile.pandas_to_dumpfile(df)#dump class processing

            else:
                #last index
                df = raw_data.loc[index:len(raw_data)+1,:]#making new dataFrame
                dump_class = dumpFile.pandas_to_dumpfile(df)#dump class processing

            #adding to dictionary
            if ids == ["TimestepDefault"]:
                #using timestep to insert
                dump_files[int(dump_class.timestep)] = dump_class
            else:
                #using custom id
                dump_files[ids[ind]] = dump_class

        return dump_files


    else:
         warnings.warn("Length of ids list is not equal to files list length")


def batch_import_files(file_paths:list,ids:list = ["TimestepDefault"]):
    """
    this opens several lammps dumps and converts it to a dictionary of
    dumpFile classes with the keys set to the timesteps

    ids:list = ["TimestepDefault"]
    ids are set to the dumpclass timestep by default however if there are duplicates
    this will override the timesteps so you can define the ids for the dictionary
    """
    if len(ids) == len(file_paths) or ids == ["TimestepDefault"]:

        dump_files = {} #dictionary of class

        for ind,file_path in enumerate(file_paths):
            #importing class
            dump_class = dumpFile.lammps_dump(file_path)

            #adding to dictionary
            if ids == ["TimestepDefault"]:
                #using timestep to insert
                dump_files[int(dump_class.timestep)] = dump_class
            else:
                #using custom id
                dump_files[ids[ind]] = dump_class

        return dump_files


    else:
         warnings.warn("Length of ids list is not equal to files list length")
