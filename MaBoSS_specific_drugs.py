#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on November 24 2017
@author: Jonas BÉAL
jonas.beal@curie.fr
"""
#%% Imports

#PC_20180802 test.txt -sy Mac -p 3 -s test -m Results/Profiles/PC_20180802_CL_mut.csv -c Results/Profiles/PC_20180802_CL_RNA_norm.csv -rb Results/Profiles/PC_20180802_CL_RNA_norm.csv -rf 100 -d Results/Profiles/PC_20180802_drug_profiles_short.csv
# Imports and tests
import sys
import os
import re
import argparse
import pandas as pd
import time
import multiprocessing as mp
import numpy
import shutil

current_wd= os.getcwd()

arg = sys.argv

parser = argparse.ArgumentParser()

#Required arguments
parser.add_argument("model", help="name of MaBoSS files, without .cfg or .bnd extension (ex: 'CancerModel')")
parser.add_argument("save_file", help="save_file is the name of the text file containing final probabilities of outputs, one simulation by line (ex: 'results.txt')")

#Optional arguments for computation parameters
parser.add_argument("-sy","--system", help="Computer OS, eiter Mac or Linux (ex: 'Mac' is you are using MacOS, 'Linux' is the default)")
parser.add_argument("-p","--num_processes", type=int, help="nb of parallel processes during simulations (ex: '3' if you want to simulate 3 profiles at the same time)")

#Optional arguments for instantation parameters

#General arguments, for all patients
parser.add_argument("-i","--inputs", help="initial probabilities of inputs, alternatively called source nodes (i.e not regulated nodes), are set to the specified value, otherwise it will be 0.5 (ex: 'Nutrients:0.3,Androgen:0.6' results in 'Nutrients = 0.3[1], 0.7[0]' and same for Androgen)")
parser.add_argument("-o","--outputs", help="outputs are marked as external nodes, whose final probabilities are saved in the result file (ex: 'Proliferation,Apoptosis')")
parser.add_argument("-s","--suffix", help="suffix is added to all intermediate and result files (ex: 'my_simulation')")

#Create a list of protected nodes, modified only based on mutations (-m) information. Therefore, initial states and rates are protected and taken from the initial .cfg file
parser.add_argument("-pn","--protected_nodes", help="list of protected nodes, modified only based on mutations (-m) information; therefore, initial states and rates are protected and taken from the initial .cfg file (ex: 'SHH,E2F4')")

#Bypass general arguments providing directly a proper CFG file with all information about inputs, outputs, define internal nodes and initial states
parser.add_argument("-cfg","--CFGbypass", type=bool, help="True if you want to ignore inputs and outputs from general arguments and extract inputs, outputs, internal nodes and initial states' information directly from the provided CFG file (ex: '-cfg True' or '-cfg 1').")

#Patient-specific arguments
parser.add_argument("-m","--mutants", help="name of the csv file containing perturbation profiles to define node mutants (also called node activity status): one profile/patient by line, with multiple genes separated by a comma. Binary 0/1 information (NA tolerated)")
parser.add_argument("-c","--init_cond", help="name of the csv file containing perturbation profiles to define node initial conditions: one profile/patient by line, with multiple genes separated by a comma. Binary 0/1 (NA tolerated) or continuous [0,1] information")
parser.add_argument("-rb","--rates_basic", help="name of the csv file containing perturbation profiles to define reaction rates based on node normalized state: one profile/patient by line, with multiple genes separated by a comma. Continuous [0,1] information")
parser.add_argument("-ra","--rates_advanced", help="name of the csv file containing perturbation profiles to define reaction rates based on activators/inhibitors nodes states: one profile/patient by line, with multiple genes separated by a comma. Continuous [0,1] information")
parser.add_argument("-rf","--rates_factor", help="multiplication factor for rates when using one of the 'rates' methods (ex: '100' in order to have rates between 1/100 and 100)")

#Patient-specific arguments
parser.add_argument("-d","--drugs", help="bla")

args = parser.parse_args()

#%% Process Arguments
print("Arguments:\n")
print(args)
print("\n")


# path to script is used to call other scripts in the same folder
pathtoscript = os.path.dirname(os.path.abspath(os.path.realpath(sys.argv[0])))
base_path = os.path.dirname(os.path.dirname(pathtoscript))+"/"
os.chdir(base_path)

#Check OS
if args.system is not None:
    system=args.system
else:
    system='Linux'

if system == 'Linux':
    sed_string="sed -i "
elif system == 'Mac':
    sed_string="sed -i '' "
else:
    print("Please use either Linux or Mac OS")
    sys.exit(1)

#Check that the model exist or translate it to maboss format if it is a .bnet file
model=args.model
path_model = base_path+"Models/"+model+"/"+model
print("Model: "+model)

if not os.path.isfile(path_model+".bnd"):
    if not os.path.isfile(path_model+".bnet"):
        print(model+".bnd and .bnet files not found")
        sys.exit(1)
    else:
        print("Convert .bnet file in .bnd file")
        sys.exit(1)
        
if not os.path.isfile(path_model+".cfg"):
    print(model+".cfg file not found")
    sys.exit(1) 
    
#Define the CFGbypass status of the simulation
if args.CFGbypass is not None:
    CFGbypass=args.CFGbypass
else:
    CFGbypass=False

#Define the protected nodes    
if args.protected_nodes is not None:
    protected_nodes = args.protected_nodes.split(",")
else:
    protected_nodes = list()

#Define all nodes, constant nodes and input nodes based on .bnd file
lines = open(path_model+".bnd").readlines()
constant_nodes = dict()
input_nodes = dict()
nodes = list()
for i in range(len(lines)):
    if re.search("^Node .*{$", lines[i]):
        node_name = re.split(" *{",lines[i])[0].split(" ")[1]        
        nodes.append(node_name)
        if re.search("logic *= *\(?[01]\)?;", lines[i+1]):
            constant_nodes[node_name] = int(lines[i+1][-3])
        if re.search("logic *= *\(?"+re.escape(node_name)+"\)?;",lines[i+1]) and (node_name not in protected_nodes):
             input_nodes[node_name] = 0.5 #Unless otherwise specified, we define a 0.5 default value for input 

#Define save_file
save_file=base_path+"Results/Simulations/"+args.save_file

#Define number of nodes
nbnodes=len(nodes)

# number of parallel processes
if args.num_processes is not None:
    nbprocesses=args.num_processes
else: nbprocesses=1
print("Number of parallel processes: "+str(nbprocesses))

#Use inputs argument to modify previous dict
if args.inputs is not None:
    inputs = args.inputs.split(",")
    for input_item in inputs:
        input_nodes[input_item.split(":")[0]] = float(input_item.split(":")[1])

#Define outputs
if os.path.isfile(path_model + ".reggraph"):
    lines = open(path_model + ".reggraph").read().splitlines()
    lines=[x.split(" ") for x in lines]
    targets = [item[0] for item in lines]
    actions = [item[1] for item in lines]
    players = [item[2] for item in lines]
    model_outputs = list(set(targets) - set(players))

if args.outputs is not None:
    outputs = args.outputs.split(",")
elif os.path.isfile(path_model + ".reggraph") and len(model_outputs)>0:
    outputs = model_outputs
else:
    outputs = ["Proliferation","Apoptosis"]

if os.path.isfile(path_model + ".reggraph"):
    print("- Inputs (model-based and user-defined): "+str(input_nodes)+"\n - Constant nodes (model-based): "+str(constant_nodes)+"\n - Simulation outputs (user-defined): "+str(outputs)+"\n - Model outputs (model-based): "+str(model_outputs))
else:
    print("- Inputs (model-based and user-defined): "+str(input_nodes)+"\n - Constant nodes (model-based): "+str(constant_nodes)+"\n - Simulation outputs (user-defined): "+str(outputs))
    
#Define patient lists
cases_common = list()

#Define mutants    
if args.mutants is not None:
    mutants = pd.read_csv(args.mutants)
    mutants.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    mutants_dict = mutants.set_index('Name').to_dict(orient="index")
    cases_mut = list(mutants_dict.keys())
    cases_common.append(cases_mut)

#Define rates
if args.rates_basic is not None:
    rates_f = float(args.rates_factor)
    rates = pd.read_csv(args.rates_basic)
    rates.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    if args.protected_nodes is not None:
        rates.drop(list(set(protected_nodes) & set(list(rates))), axis=1, inplace=True)
    rates_dict = rates.set_index('Name').to_dict(orient="index")
    cases_rates = list(rates_dict.keys())
    cases_common.append(cases_rates)
    
#Define rates_advanced and factor
if args.rates_advanced is not None:
    
    #Import data to define rates
    rates_a = pd.read_csv(args.rates_advanced)
    rates_a.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    if args.protected_nodes is not None:
        rates_a.drop(list(set(protected_nodes) & set(list(rates_a))), axis=1, inplace=True)
    rates_a_dict = rates_a.set_index('Name').to_dict(orient="index")
    cases_rates_a = list(rates_a_dict.keys())
    cases_common.append(cases_rates_a)
    rates_f = float(args.rates_factor)
    
    rates_nodes=list(rates_a)
    rates_nodes.remove('Name')
    
    #Activator/Inhibitor dictionnaries
    model_activators = {}
    for node in set(targets):
        model_activators[node] = [x for x,y,z in zip(players,actions,targets) if (x in rates_nodes) and y=='->' and z==node]
    model_inhibitors = {}
    for node in set(targets):
        model_inhibitors[node] = [x for x,y,z in zip(players,actions,targets) if (x in rates_nodes) and y=='-|' and z==node]
    
    model_activators={k: v for k, v in model_activators.items() if v}
    model_inhibitors={k: v for k, v in model_inhibitors.items() if v}
    
#Define init_cond profile
if args.init_cond is not None:
    init_cond = pd.read_csv(args.init_cond)
    init_cond.rename(columns={'Unnamed: 0':'Name'}, inplace=True)
    if args.protected_nodes is not None:
        init_cond.drop(list(set(protected_nodes) & set(list(init_cond))), axis=1, inplace=True)
    init_cond_dict = init_cond.set_index('Name').to_dict(orient="index")
    cases_exp = list(init_cond_dict.keys())
    cases_common.append(cases_exp)

#Define drugs profile
drugs_dict = {}
drugs_dict['Zero'] = {}
if args.drugs is not None:
    drugs = pd.read_csv(args.drugs)
    for i in drugs['DRUG_ID'].unique():
        drugs_dict[i] = {drugs['TARGETED_NODE'][j]: drugs['STRENGTH'][j] for j in drugs[drugs['DRUG_ID']==i].index}        
    nbnodes+=len(drugs_dict)-1
    
#MaBoSS exec
if nbnodes<=64:
    maboss_exec = base_path+"MaBoSS/"+system+"/MaBoSS"
elif nbnodes<=150:
    maboss_exec = base_path+"MaBoSS/"+system+"/MaBoSS_150n"
else:
    print("Your model has more than 150 nodes, please recompile MaBoSS with the number of nodes of your model. See MaBoSS compiling help: http://maboss.curie.fr/")
if not os.path.isfile(maboss_exec):
    print("Relevant MaBoSS executable is not available")
    sys.exit(1)
print("MaBoSS executable: "+maboss_exec)
os.system("chmod +x "+maboss_exec)

#Define patient lists
if not all(v is None for v in [args.mutants, args.rates_basic, args.rates_advanced, args.init_cond]):
    cases_common = list(set.intersection(*map(set, cases_common)))

#Define suffix
if args.suffix is not None:
    suffix = args.suffix
else:
    suffix = "classic"
print("Suffix: "+suffix)
    

#%% Start modifying files



#Prepare simulation copying model files for subsequent modifications
#os.system("cp "+path_model+".bnd "+path_model+"_"+suffix+".bnd")
shutil.copy2(path_model+".bnd", path_model+"_"+suffix+".bnd")
#os.system("cp "+path_model+".cfg "+path_model+"_"+suffix+".cfg")
shutil.copy2(path_model+".cfg", path_model+"_"+suffix+".cfg")

model=model+"_"+suffix
path_model=path_model+"_"+suffix

#Define outputs as external nodes
if not CFGbypass:
    for node in nodes:
        os.system(sed_string+"'/^\[*"+node+"\]*\.is_internal *= */d' "+path_model+".cfg")
        if node in outputs:
            os.system("echo '"+node+".is_internal = FALSE;' >> "+path_model+".cfg")
        else:
            os.system("echo '"+node+".is_internal = TRUE;' >> "+path_model+".cfg")

#Define proper initial conditions for inputs and constant nodes (implicit inputs)
if not CFGbypass:
    for input_item, input_value in dict(input_nodes, **constant_nodes).items():
        os.system(sed_string+"'/^\[*"+input_item+"\]*\.istate *=/d' "+path_model+".cfg")
        os.system("echo '["+input_item+"].istate = "+str(input_value)+"[1], "+str(1-input_value)+"[0];' >> "+path_model+".cfg")

#%%Define function used to perform the simulation itself and process the output
def perform_MaBoSS_simulation(n_profile, fname, profile_name, list_drugs):
    if profile_name is not "WT":
        path_fname=path_model+"_"+str(n_profile)
    else:
        path_fname=path_model
    
    for drug_id in list_drugs:
        path_fname_drug=path_fname+"_drug"+str(drug_id)
        fname_drug=fname+"_drug"+str(drug_id)
        string_opt=path_fname_drug+".bnd " + path_fname_drug + ".cfg -mb " + maboss_exec
        os.system("chmod +x "+path_fname_drug+".bnd")
        os.system("chmod +x "+path_fname_drug+".cfg")
        os.chdir(os.path.dirname(path_model))
        os.system("perl "+base_path+"Scripts/Simulations/MBSS_FormatTable.pl " + string_opt)
        os.chdir(current_wd)
        
        # store column names (states) and last line with final state distribution in a temp file
        os.system("head -n 1 "+path_fname_drug+"/"+fname_drug+"_probtraj_table.csv > "+path_fname_drug+"/"+fname_drug+"_lastprob_table.csv")
        os.system("tail -n 1 "+path_fname_drug+"/"+fname_drug+"_probtraj_table.csv >> "+path_fname_drug+"/"+fname_drug+"_lastprob_table.csv")
        # extract the output probas from final state distribution in fname_drug+"/"+fname_drug+"_lastprob.csv
        #os.system("Rscript "+base_path+"Scripts/Simulations/extract_output_probas.R -i "+path_fname_drug+"/"+fname_drug+"_lastprob_table.csv -o "+path_fname_drug+"/"+fname_drug+"_lastprob -n "+",".join(outputs))
        os.system("/bioinfo/local/build/Centos/R/R-3.4.0/bin/Rscript "+pathtoscript+"/extract_output_probas.R -i "+path_fname_drug+"/"+fname_drug+"_lastprob_table.csv -o "+path_fname_drug+"/"+fname_drug+"_lastprob -n "+",".join(outputs))
    	# add the number of the simulation to the line finally saved
        os.system("echo "+str(n_profile)+"\t"+profile_name+"\t"+str(drug_id)+" $(tail -n 1 "+path_fname_drug+"/"+fname_drug+"_lastprob.csv) >> "+save_file)
        #Remove temporary files
        os.system("rm "+path_fname_drug+"/*")
        os.system("rmdir "+path_fname_drug)
        os.system("rm "+path_fname_drug+".bnd")
        os.system("rm "+path_fname_drug+".cfg")
    if profile_name is not "WT":
        previous.append(n_profile)

#%%Launch simulations depending on the case
os.system("echo n_profile\tPatient_ID\tDrug_ID\tTime\t"+'\t'.join(outputs)+"\tTH > "+save_file)
if all(v is None for v in [args.mutants, args.rates_basic, args.rates_advanced, args.init_cond]):
    n_profile=0
    fname=model
    perform_MaBoSS_simulation(0,fname, "WT")
    
else:
    manager=mp.Manager()
    previous=manager.list()
    processes = list()
    for i in range(len(cases_common)):
        patient_id=cases_common[i]
        fname=model+"_"+str(i)
        path_fname=path_model+"_"+str(i)
        #os.system("cp "+path_model+".bnd "+path_fname+".bnd")
        shutil.copy2(path_model+".bnd", path_fname+".bnd")
        #os.system("cp "+path_model+".cfg "+path_fname+".cfg")
        shutil.copy2(path_model+".cfg", path_fname+".cfg")
    
        # set init_cond profiles        
        if args.init_cond is not None:
            patient_dict = init_cond_dict[patient_id]
            patient_dict_red = { k:v for k, v in patient_dict.items() if not numpy.isnan(v) }
            for node, value in patient_dict_red.items():
                value_red = round(value,5)
                os.system(sed_string+"'/^\[*"+node+"\]*\.istate/d' "+path_fname+"'.cfg'")
                os.system("echo '["+node+"].istate = "+str(value_red)+"[1], "+str(1-value_red)+"[0];' >> "+path_fname+".cfg")
        
        # set rates profiles        
        if args.rates_basic is not None:
            rates_list=rates_dict[patient_id]
            rates_list_red={ k:v for k, v in rates_list.items() if not numpy.isnan(v) }
            for node, value in rates_list_red.items():
                if not numpy.isnan(value):
                    up_value = round(rates_f**(2*value-1),5)
                    down_value = round(1/up_value,5)
                    original_up = float(os.popen("grep -E '^\$u_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    original_down = float(os.popen("grep -E '^\$d_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    os.system(sed_string+"'s/u_"+node+" *= *[0-9]*\.*[0-9]*;/u_"+node+"="+str(up_value*original_up)+";/g' "+path_fname+"'.cfg'")
                    os.system(sed_string+"'s/d_"+node+" *= *[0-9]*\.*[0-9]*;/d_"+node+"="+str(down_value*original_down)+";/g' "+path_fname+"'.cfg'")
                    value_red = round(value,5)
                    os.system(sed_string+"'/^\[*"+node+"\]*\.istate/d' "+path_fname+"'.cfg'")
                    os.system("echo '["+node+"].istate = "+str(value_red)+"[1], "+str(1-value_red)+"[0];' >> "+path_fname+".cfg")
                    
        # set rates_advanced profiles         
        if args.rates_advanced is not None:
            rates_dict_patient=rates_a_dict[patient_id]
            
            for node, value in model_activators.items():
                acti_value = numpy.nanmean([rates_dict_patient[node_name] for node_name in value])
                if not numpy.isnan(acti_value):
                    original_up = float(os.popen("grep -E '^\$u_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    rate_value = round(rates_f**(2*acti_value-1),5)
                    os.system(sed_string+"'s/u_"+node+" *= *[0-9]*\.*[0-9]*;/u_"+node+"="+str(rate_value*original_up)+";/g' "+path_fname+"'.cfg'")
            
            for node, value in model_inhibitors.items():
                inhi_value = numpy.nanmean([rates_dict_patient[node_name] for node_name in value])
                if not numpy.isnan(inhi_value):
                    original_down = float(os.popen("grep -E '^\$d_"+node+" *= *' "+path_fname+".cfg | cut -d'=' -f2 | cut -d';' -f1| tr -d '\n'").read())
                    rate_value = round(rates_f**(2*inhi_value-1),5)
                    os.system(sed_string+"'s/d_"+node+" *= *[0-9]*\.*[0-9]*;/d_"+node+"="+str(rate_value*original_down)+";/g' "+path_fname+"'.cfg'")
            
            inputs_to_specify = list(set(input_nodes.keys()) & set(rates_nodes))
            for node in inputs_to_specify:
                value = round(rates_dict_patient[node],5)
                os.system(sed_string+"'/^\[*"+node+"\]*\.istate/d' "+path_fname+"'.cfg'")
                os.system("echo '["+node+"].istate = "+str(value)+"[1], "+str(1-value)+"[0];' >> "+path_fname+".cfg")

        # set mutants profiles 
        if args.mutants is not None:
            mutants_list=mutants_dict[patient_id]
            mutants_list_red={ k:v for k, v in mutants_list.items() if not numpy.isnan(v) }
            for node, value in mutants_list_red.items():
                os.system(sed_string+"'/^\[*"+node+"\]*\.istate/d' "+path_fname+"'.cfg'")
                if value==0:
                    os.system(sed_string+"'s/u_"+node+" *= *[0-9]*\.*[0-9]*;/u_"+node+"=0;/g' "+path_fname+"'.cfg'")
                    os.system("echo '["+node+"].istate=0[1], 1[0];' >> "+path_fname+".cfg")
                elif value==1:
                    os.system(sed_string+"'s/d_"+node+" *= *[0-9]*\.*[0-9]*;/d_"+node+"=0;/g' "+path_fname+"'.cfg'")
                    os.system("echo '["+node+"].istate=1[1], 0[1];' >> "+path_fname+".cfg")
        
        # implement drug effects
        if args.drugs is not None:
            file_content = open(str(path_fname+".bnd")).read().replace("\n","")
            bnd_content = re.sub(r"(Node +[^{]+)",r"\n\1", file_content)
            cfg_content = open(str(path_fname+".cfg")).read()
            
            #Modify the BND file adding inhibitor nodes and their actions
            for target in set(drugs.TARGETED_NODE):
                bnd_content = re.sub(r"(Node +"+target+" *{ *logic *= *)([^;]+)", r"\1(\2) & !anti_"+target,bnd_content)
                bnd_content += "\nNode anti_{}".format(target) + " { " + "logic = (anti_{});".format(target) + " rate_up = @logic ? 0 : 0; rate_down = @logic ? 0 : 0;}"
                cfg_content += "\nanti_"+target+".is_internal = TRUE;\n[anti_"+target+"].istate = 0[1], 1[0];"
            bnd_content = bnd_content.replace("{", "{\n").replace("}", "}\n").replace(";", ";\n")
            
            #Generate drug specific CFG files and their BND counterparts
            for drug_id, drug_indiv_dict in drugs_dict.items():
                drug_cfg_content=cfg_content
                for target, strength in drug_indiv_dict.items():
                    value=round(strength,3)
                    drug_cfg_content = re.sub(r"(\[anti_"+target+"\].istate *= *)[^;]+", r"\1 "+str(value)+" [1], "+str(1-value)+"[0]", drug_cfg_content)
                open(path_fname+"_drug"+str(drug_id)+".cfg","w").write(drug_cfg_content)
                open(path_fname+"_drug"+str(drug_id)+".bnd","w").write(bnd_content)
        else:
            os.rename(path_fname+".bnd", path_fname+"_drugZero.bnd")
            os.rename(path_fname+".cfg", path_fname+"_drugZero.cfg")
        
        while len(previous)<i-(nbprocesses-1):
            time.sleep(1)
        print(str(i)+": "+patient_id)
        p = mp.Process(target = perform_MaBoSS_simulation, args=(i,fname,patient_id, list(drugs_dict)))
        p.start()
        processes.append(p)
    
    for process in processes:
        process.join()
        
os.system("rm "+path_model+".bnd")
os.system("rm "+path_model+".cfg")
