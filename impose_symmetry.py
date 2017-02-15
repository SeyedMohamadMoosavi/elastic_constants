#! /usr/bin/env python
import numpy as np
import spglib
import math
import sys
from symm_functions import *
################### INPUT  ##########################################
str_name=sys.argv[1]
precision=1
cell,types=read_cif(str_name)
n_at_o=len(cell[2])
######################################################################
spcg=spglib.get_spacegroup(cell, symprec=precision)
print(spcg)
lattice_new, scaled_positions_new, numbers_new = spglib.standardize_cell(cell, to_primitive=False, no_idealize=False, symprec=precision,angle_tolerance=precision)
# lattice_new, scaled_positions_new, numbers_new = spglib.find_primitive(cell, symprec=precision,angle_tolerance=1)
cell_new=(lattice_new,scaled_positions_new,numbers_new)
n_at_n=len(numbers_new)
print_cif('primitive_'+str_name,cell_new,types)
f_info=open("INFO",'a')
f_info.write("%s %10i %10i %10s\n"%(str_name,n_at_o, n_at_n ,spcg))
