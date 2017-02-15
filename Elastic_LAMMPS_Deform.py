#! /usr/bin/env python
import math
import argparse

######
# input file name
######
parser = argparse.ArgumentParser("using cell parameters from a cif file to deform a structure in lammps ")
parser.add_argument("--filecif","-fc", help ="name of the cif file")
parser.add_argument("--filein","-fi", help ="name of the input file")
args = parser.parse_args()


bohr2ang=0.529177249
DEG2RAD=math.pi/180.
#############

cif_name=args.filecif
cond=True
cond2=False
f_cif=open(cif_name,'r')
atom_props_count=0
atoms=[]
for line in f_cif.readlines(): 
    line_stripped=line.strip()
    if (not line) or line_stripped.startswith("#"):
        continue
    line_splitted=line.split()
    if line_stripped.startswith("_cell_length_a"):
        cell_a=float(line_splitted[1])
    elif line_stripped.startswith("_cell_length_b"):
        cell_b=float(line_splitted[1])
    elif line_stripped.startswith("_cell_length_c"):
        cell_c=float(line_splitted[1])
    elif line_stripped.startswith("_cell_angle_alpha"):
        cell_alpha=float(line_splitted[1])
    elif line_stripped.startswith("_cell_angle_beta"):
        cell_beta=float(line_splitted[1])
    elif line_stripped.startswith("_cell_angle_gamma"):
        cell_gamma=float(line_splitted[1])
    if cond2==True and line_stripped.startswith("loop_"):
        break
    else:
        if line_stripped.startswith("_atom") :
            atom_props_count+=1
            if line_stripped=="_atom_site_type_symbol":
                type_index=atom_props_count-1
            elif line_stripped=="_atom_site_fract_x":
                fracx_index=atom_props_count-1 
            elif line_stripped=="_atom_site_fract_y":
                fracy_index=atom_props_count-1
            elif line_stripped=="_atom_site_fract_z":
                fracz_index=atom_props_count-1

            cond2=True
        elif cond2==True:
            if len(line_splitted)==atom_props_count:
                atoms.append(line)
f_cif.close()
lx = cell_a
xy = cell_b*math.cos(cell_gamma*DEG2RAD)
xz = cell_c*math.cos(cell_beta*DEG2RAD)
ly = math.sqrt(cell_b**2 - xy**2)
yz = (cell_b*cell_c*math.cos(cell_alpha*DEG2RAD) - xy*xz)/ly
lz = math.sqrt(cell_c**2 - xz**2 - yz**2)

input_name=args.filein
f_input=open(input_name,'r')
data=f_input.readlines()

for i,line in enumerate(data):
    if ("read_dump" and "relaxed_structure.dump") in line:
        index_read_data=i
        break

strin="change_box all "
strin+=" x final %10.8f %10.8f "%(0.0,lx)  
strin+=" y final %10.8f %10.8f "%(0.0,ly)
strin+=" z final %10.8f %10.8f "%(0.0,lz)
strin+=" xy final %10.8f"%(xy) 
strin+=" xz final %10.8f"%(xz)
strin+=" yz final %10.8f"%(yz)
strin+=" remap units box"
strin+="\n"
data.insert(index_read_data+1,strin)
f_input.close()
f_out=open(input_name,'w')
for line in data:
    f_out.write(line)
f_out.close()

