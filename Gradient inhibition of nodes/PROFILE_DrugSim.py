#!/usr/bin/env python3
# coding: utf-8
#
# MaBoSS (Markov Boolean Stochastic Simulator)
# Copyright (C) 2011-2017 Institut Curie, 26 rue d'Ulm, Paris, France
#   
# MaBoSS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#   
# MaBoSS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#   
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA 
#
"""
Created on May 2018
Last edit on June 2021
@author: Arnau Montagud, arnau.montagud@gmail.com
Built on code by Barthélémy Caron and Jonas Béal
"""
#%% Imports
import os, sys
import subprocess
import time
import numpy as np
from os import listdir
import shutil
import re
import argparse

current_wd= os.getcwd()
arg = sys.argv
parser = argparse.ArgumentParser()

#Required arguments
parser.add_argument("model", help="Name of MaBoSS files, without .cfg or .bnd extension (ex: 'CancerModel').")
parser.add_argument("-d","--drugs", help="Names of nodes affected by drug, commma separated (ex: 'PI3K, MEK').")
parser.add_argument("-c","--drug_conc", help="Levels of drug inhibition, between 0 and 1 (ex: '0, 0.2, 0.4, 0.6, 0.8, 1').")
# example: LNCAP_AR_EGF -d "PI3K, MEK" -c "0, 0.2, 0.4, 0.6, 0.8, 1"

args = parser.parse_args()

#%% Process Arguments
print("Arguments:\n")
print(args)
print("\n")

#Check that the model exist or translate it to maboss format if it is a .bnet file
model=args.model
print("Model: "+model)

if not os.path.isfile(model+".bnd"):
    print(model+".bnd file not found")
    sys.exit(1)
    
if not os.path.isfile(model+".cfg"):
    print(model+".cfg file not found")
    sys.exit(1)

drugs=args.drugs
node_list = str(drugs.replace(" ","").replace(",",", ")).split(", ")
print("Nodes to be inhibited: "+str(node_list).replace("[","").replace("]","").replace("'", ""))

text = '2,4,6,8|10,12,14,16|18,20,22,24'
my_data = [x.split(',') for x in text.split('|')]
my_data
[['2', '4', '6', '8'], ['10', '12', '14', '16'], ['18', '20', '22', '24']]

BNDnodes1 = [re.findall(r'Node.*',line) for line in open(str(model+".bnd"))]
BNDnodes2 = list(filter(None, BNDnodes1))
BNDnodes = [[w.replace("Node ", "").replace(" {", "") for w in line] for line in BNDnodes2]
node_list1 = [x.split(", ") for x in node_list]
drugs_notin_BND = [item for item in node_list1 if item not in BNDnodes]
if len(drugs_notin_BND) > 0:
    print("Nodes NOT in BND file: "+str(drugs_notin_BND).replace("[","").replace("]","").replace("'", ""))

drug_conc=args.drug_conc
istate_value_list = [float(i) for i in drug_conc.replace(" ","").replace(",",", ").split(", ")]
istate_value_name_list = [str(i).replace(".", "_") for i in istate_value_list]
print("Inhibited nodes' levels: "+str(istate_value_list).replace("[","").replace("]",""))

output= str("drugsim_"+model+"_"+str(node_list).replace("[","").replace("]","").replace("'", "").replace(", ", "_"))

if os.path.exists(output):
    shutil.rmtree(output)
os.mkdir(output)

#%% Generate new BND file containing the inhibitor nodes and altering the logical rules corresponding to the inhibitions
print("Modifying BND file: adding inhibitor nodes and inhibitions")
file_content = open(str(model+".bnd")).read()
file_content = file_content.replace(" {", "{")
file_content = file_content.replace("\n{", "{")
file_content = file_content.replace("{", "\n{")

node_list_from_file = file_content.split("Node")[1:]
for item_index, item in enumerate(node_list_from_file):
    node_list_from_file[item_index] = "Node" + item

new_node_list = str()
for node in node_list:
    mystring = "Node {}\n".format(node)
    item_count = 0

    for item_index, item in enumerate(node_list_from_file):
        if mystring in (node_list_from_file[item_index]):
            item_count += 1
            for line_index, line in enumerate(node_list_from_file[item_index].split("\n")):
                if " logic" in line:
                    new_node_list += line.split(";")[0].split(" = ")[0] + " = (" + line.split(";")[0].split(" = ")[1] + ") AND NOT anti_{};\n".format(node)
                else:
                    new_node_list += line + "\n"
            node_list_from_file.pop(item_index)
    if item_count == 0:
        print("{} not found".format(node))

for last_item in node_list_from_file:
    new_node_list += last_item

for inhibitor_to_add in node_list:
    new_node_list += "Node anti_{}".format(inhibitor_to_add) + " \n{\n" + "\tlogic = (anti_{});\n".format(inhibitor_to_add) + "\trate_up = @logic ? 0 : 0;\n\trate_down = @logic ? 0 : 0;\n}\n\n\n"
    
name_bnd = "{}_{}_new.bnd".format(model,"".join(node_list))
fw1 = open(name_bnd, "w")
fw1.write(new_node_list)
fw1.close()

#%%Generate new CFG file setting the initial parameters of the new nodes
print("Modifying CFG file: istate formatting")
file_content = open(str(model+".cfg")).read()

file_content = file_content.replace("[1]", " [1]").replace("[0]", " [0]").replace("] ,", "],").replace("=", " = ").replace("  ", " ")

new_cfg = str()
for line_index, line in enumerate(file_content.split("\n")):
    if ".istate = TRUE" in line:
        line = "[" + line.replace(".istate = TRUE", "].istate = 1 [1], 0 [0]")
    elif ".istate = FALSE" in line:
        line = "[" + line.replace(".istate = FALSE", "].istate = 0 [1], 1 [0]")
    elif "].istate = 0 [0], 1 [1]" in line:
        line = line.replace("].istate = 0 [0], 1 [1]", "].istate = 1 [1], 0 [0]")
    elif "].istate = 1 [0], 0 [1]" in line:
        line = line.replace("].istate = 1 [0], 0 [1]", "].istate = 0 [1], 1 [0]")
    else:
        line = line
    new_cfg += line + "\n"  

b1 = [re.findall(r'.*istate.*',line) for line in open(str(model+".cfg"))]
b2 = list(filter(None, b1))
CFGnodes = [[re.sub(r"].istate.*","",w) for w in line] for line in b2]
CFGnodes = [[re.sub(r"\[","",w) for w in line] for line in CFGnodes]

a4 = [item for item in BNDnodes if item not in CFGnodes]
if len(a4) > 0:
    a5 = str([re.sub(r"$",".istate = 0.5 [0] , 0.5 [1];",str(line)) for line in a4]).replace("'","").replace('"[',"").replace('"]',"").replace(";",";\n").replace('", ','[')
    new_cfg += a5

name_cfg = "{}_{}_new.cfg".format(model,"".join(node_list))
fw2 = open(name_cfg, "w")
fw2.write(new_cfg)
fw2.close()

#%% Create BND and CFG files with the required inhibition levels when at least two drugs are specified
if len(node_list) != 1:
    for node_index_1, node_1 in enumerate(node_list):
        for node_index_2, node_2 in enumerate(node_list[1+node_index_1:]):
            print("creating {} {} combinations".format(node_1, node_2))
            for istate_index_1, value_1 in enumerate(istate_value_list):
                for istate_index_2, value_2 in enumerate(istate_value_list):
                    value_1_name = istate_value_name_list[istate_index_1]
                    value_2_name = istate_value_name_list[istate_index_2]
                    """GEN BND FILE"""
                    options_copy = [name_bnd, 
                                    str(output+"/{}_{}_{}_{}_{}.bnd").format(
                                            model, 
                                            node_1, 
                                            istate_value_name_list[istate_index_1], 
                                            node_2, 
                                            istate_value_name_list[istate_index_2])]
                    shutil.copy(options_copy[0], options_copy[1])
                    """GEN CFG FILE"""
                    f = open(name_cfg,"r")
                    lines = f.readlines()
                    f.close()
                    save_name = str(output+"/{}_{}_{}_{}_{}.cfg").format(
                            model, 
                            node_1, 
                            value_1_name, 
                            node_2, 
                            value_2_name)
                    nf = open(save_name, "w")

                    for line in lines:
                        if "[{}].istate = ".format(node_1) in line:
                            new_line = line
                            new_line = "[{}].istate = {} [1], {} [0];\n[anti_{}].istate = {} [1], {} [0];\n".format(node_1, np.round(1-value_1, 1), value_1,node_1, value_1, np.round(1-value_1, 1))
                            nf.write(str(new_line))
                        elif "[{}].istate = ".format(node_2) in line:
                            new_line = line
                            new_line = "[{}].istate = {} [1], {} [0];\n[anti_{}].istate = {} [1], {} [0];\n".format(node_2, np.round(1-value_2, 1), value_2,node_2, value_2, np.round(1-value_2, 1))
                            nf.write(str(new_line))
                        else:
                            nf.write(str(line))
                    nf.close()


#Create BND and CFG files with the required inhibition levels when only one drug has been added
else:
    for node_index_1, node_1 in enumerate(node_list):
        print("creating {} gradual inhibition".format(node_1))
        for istate_index_1, value_1 in enumerate(istate_value_list):
            value_1_name = istate_value_name_list[istate_index_1]
            """GEN BND FILE"""
            options_copy = [name_bnd, 
                                    str(output+"/{}_{}_{}.bnd").format(model, node_1, istate_value_name_list[istate_index_1])]
            proc2 = subprocess.Popen(["cp"] + options_copy)
            time.sleep(1.0)
            """GEN CFG FILE"""
            f = open(name_cfg,"r")
            lines = f.readlines()
            f.close()
            save_name = str(output+"/{}_{}_{}.cfg").format(
                    model, 
                    node_1, 
                    value_1_name)
            nf = open(save_name, "w")
            for line in lines:
                if "[{}].istate = ".format(node_1) in line:
                    new_line = line
                    new_line = "[{}].istate = {} [1], {} [0];\n[anti_{}].istate = {} [1], {} [0];\n".format(node_1, np.round(1-value_1, 1), value_1,node_1, value_1, np.round(1-value_1, 1))
                    nf.write(str(new_line))    
                else:
                    nf.write(str(line))
            nf.close()

file_list = []
count = 0
for file_index, file in enumerate(listdir(output)):
    if ".bnd" in file:
        file_list.append(file)
number_of_simulations = len(file_list)

#%% output a runnning script

if len(file_list) != 0:
    fw1 = open("./run_"+output+".sh", "w")
    fw1.write("#!/bin/bash\n" )
    fw1.write(str("cd "+output+"\n"))
    for item_list, item in enumerate(file_list):
        item = item.strip(".bnd")
        fw1.write("echo {}/{} mutant: {}\n".format(item_list+1,len(file_list),item))
        fw1.write("../MBSS_FormatTable.pl " + "{}.bnd ".format(item) + "{}.cfg ".format(item) + "\n")
        fw1.write("rm -f *_fp.csv *_probtraj.csv *_run.txt *_statdist.csv\n")
    fw1.close()

os.chmod("./run_"+output+".sh", 0o755)
os.system("./run_"+output+".sh")
os.system("rm -f ./run_"+output+".sh ./"+name_bnd+" ./"+name_cfg)
