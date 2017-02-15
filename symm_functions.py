#! /usr/bin/env python
import numpy as np
import spglib
import math
deg2rad = math.pi/180.0
########## Functions ################# 
def print_cif(file_name,cell,atom_labels):
        lattice_vector=cell[0]
        atomic_positions=cell[1]
        atom_types=cell[2]
        fo = open(file_name , 'w')
        ## writing the general P1 symmetry ##
        fo.write("data_%s\n"%file_name)
        fo.write("_audit_creation_method           Elastic_mohamad\n")
        fo.write("_symmetry_space_group_name_H-M    P1\n")
        fo.write("_symmetry_Int_Tables_number       1\n")
        fo.write("_symmetry_cell_setting            triclinic\n")
        fo.write("loop_\n")
        fo.write("_symmetry_equiv_pos_as_xyz\n")
        fo.write("\'x, y, z\'\n")
        cell_a=math.sqrt(lattice_vector[0][0]**2+lattice_vector[0][1]**2+lattice_vector[0][2]**2)
        cell_b=math.sqrt(lattice_vector[1][0]**2+lattice_vector[1][1]**2+lattice_vector[1][2]**2)
        cell_c=math.sqrt(lattice_vector[2][0]**2+lattice_vector[2][1]**2+lattice_vector[2][2]**2)
        cell_gamma=math.acos(lattice_vector[1][0]/cell_b)
        cell_beta =math.acos(lattice_vector[2][0]/cell_c)
        cell_alpha=math.acos(lattice_vector[2][1]/cell_c*math.sin(cell_gamma)+math.cos(cell_beta)*math.cos(cell_gamma))
        cell_gamma=np.degrees( cell_gamma)
        cell_beta =np.degrees( cell_beta )
        cell_alpha=np.degrees( cell_alpha)
        fo.write("_cell_length_a %12.08f\n"%(cell_a))
        fo.write("_cell_length_b %12.08f\n"%(cell_b))
        fo.write("_cell_length_c %12.08f\n"%(cell_c))
        fo.write("_cell_angle_alpha %12.08f\n"%(cell_alpha))
        fo.write("_cell_angle_beta  %12.08f\n"%(cell_beta ))
        fo.write("_cell_angle_gamma %12.08f\n"%(cell_gamma))
        fo.write("loop_\n")
        fo.write("_atom_site_type_symbol\n")
        fo.write("_atom_site_fract_x\n")
        fo.write("_atom_site_fract_y\n")
        fo.write("_atom_site_fract_z\n")
        fo.write("_atom_site_label\n")
        i=1
        for at_type , at_position in zip(atom_types,atomic_positions):
            at_loc=""
            for loc in at_position:
                at_loc+="%10.4f "%loc
            fo.write("%-10s %s %10s%i\n"%(atom_labels[at_type-1],at_loc,atom_labels[at_type-1],i))
            i+=1

def read_cif(name):
    fi = open(name, 'r')
    EIF = fi.readlines()
    fi.close()

    cond=True
    cond2=False
    atom_props_count=0
    atoms=[]
    counter=0
    cell_parameter_boundary=[0.0,0.0]
    for line in EIF:
        line_stripped=line.strip()
        if (not line) or line_stripped.startswith("#"):
            continue
        line_splitted=line.split()
        if line_stripped.startswith("_cell_length_a"):
            cell_a=line_splitted[1]
            cell_parameter_boundary[0]=counter+1
        elif line_stripped.startswith("_cell_length_b"):
            cell_b=line_splitted[1]
        elif line_stripped.startswith("_cell_length_c"):
            cell_c=line_splitted[1]
        elif line_stripped.startswith("_cell_angle_alpha"):
            cell_alpha=line_splitted[1]
        elif line_stripped.startswith("_cell_angle_beta"):
            cell_beta=line_splitted[1]
        elif line_stripped.startswith("_cell_angle_gamma"):
            cell_gamma=line_splitted[1]
            cell_parameter_boundary[1]=counter+1
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
    
        counter+=1
    
    positions=[]
    numbers=[]
    types=[]
    for at in atoms:
        ln=at.strip().split()
        positions.append([float(ln[fracx_index]),float(ln[fracy_index]),float(ln[fracz_index])])
        types.append(ln[type_index])
        
    unique_types=list(set(types))
    
    for label in types:
        numbers.append(int(unique_types.index(label)+1))
    cpar=[cell_a,cell_b,cell_c]
    cang=[cell_alpha,cell_beta,cell_gamma]
    cpar=[float(i) for i in cpar]
    cang=[float(i)*deg2rad for i in cang]
    v=math.sqrt(1-math.cos(cang[0])**2-math.cos(cang[1])**2-math.cos(cang[2])**2+2*math.cos(cang[0])*math.cos(cang[1])*math.cos(cang[2]))
    CP=[]
    CP.append([cpar[0],0,0])
    CP.append([cpar[1]*math.cos(cang[2]),cpar[1]*math.sin(cang[2]),0])
    CP.append([cpar[2]*math.cos(cang[1]),cpar[2]*(math.cos(cang[0])-math.cos(cang[1])*math.cos(cang[2]))/(math.sin(cang[2])),cpar[2]*v/math.sin(cang[2])])
    CP=np.array(CP)
    cell=(CP,positions,numbers)
    return cell, unique_types
