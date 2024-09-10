import os
import sys
from typing import Optional, Dict, List, Tuple, Union, Any
from recur.utils import util
import shutil
from recur import __version__, recur_iqtree_nthreads, iqtree_nthreads, recur_nthreads
from recur import helpinfo
import fnmatch
import dendropy
import multiprocessing as mp     


class Options(object):
    
    def __init__(self):
        self.nthreads = recur_nthreads
        self.recur_nthreads, self.iqtree_nthreads = recur_iqtree_nthreads, iqtree_nthreads
        self.tree_program = "iqtree2"
        self.sequence_type = "CODON1"
        self.evolution_model = "TEST"
        self.qStartFromMSA = False
        self.iqtree_version = "local"
        self.gene = None
        self.alnpre = None
        self.name = ""   # name to identify this set of results
        self.extended_filename = False
        self.output_prefix = None
        self.gene_tree = None
        self.root_node = None
        self.nalign = 1000
        self.mcs_batch_size = None
        self.outgroups = None 
        self.bootstrap = 1000
        self.sh_alrt = 1000
        self.usr_node_aln = None 
        self.usr_state = None 
        self.usr_tree = None
        self.usr_iqtree = None
        self.usr_mcs_treeDir = None 
        self.usr_mcs_faDir = None 
        self.usr_mcs_alnDir = None
        self.keepprev = False
        self.recDir = None
        self.restart_from = 0
        self.fd_limit = None
        self.seed = 8
        self.mcs_seed = self.seed
        self.fix_branch_length = False
        self.update_cycle = None
        self.system_info = False
        self.override = True

    def what(self) -> None:
        for k, v in self.__dict__.items():
            if v == True:
                print(k)

def GetDirectoryArgument(arg: str) -> str:
    directory = os.path.abspath(arg)
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        util.Fail()
    if not os.path.isfile(directory) and directory[-1] != os.sep: 
        directory += os.sep
    return directory

def GetFileArgument(arg: str) -> str:
    file_path = os.path.abspath(arg)
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        print("Directory points to the file doesn't exist: %s" % directory)
        util.Fail()
    if not os.path.isfile(file_path):
        print("Specified file doesn't exist: %s" % file_path)
        util.Fail()
    return file_path

def validate_newick_tree(tree_file_path):
    try:
        tree = dendropy.Tree.get(path=tree_file_path, schema="newick")
        return True
    except Exception as e:
        print(f"ERROR reading tree file: {e}")
        print("Unrecognized tree file format detected. The file must be in Newick format.")
        return False

def ProcessArgs(args: List[Any]) -> Tuple[Options, str, Optional[str], Optional[str]]:
    
    options = Options()
    alnDir = ""
    alnPath = None 
    resultsDir_nonDefault = None
    usr_iqtree_nthread = False 
    usr_recur_nthread = False
    
    while len(args) > 0:
        arg = args.pop(0)

        if arg == "-f" or arg == "--faln":
            if options.qStartFromMSA:
                print("Repeated argument: -f/--faln\n")
                util.Fail()

            options.qStartFromMSA = True

            if len(args) == 0:
                print("Missing option for command line argument %s" % arg)
                util.Fail()

            isfile = False
            isdir = False
            arg = args.pop(0)

            if os.path.isfile(arg):
                isfile = True
                alnPath = GetFileArgument(arg)
                alnDir = os.path.dirname(alnPath) 
                if "." in alnPath:
                    options.gene = os.path.basename(alnPath).split(".", 1)[0]
                else:
                    options.gene = os.path.basename(alnPath)
            elif os.path.isdir(arg):
                isdir = True
                alnDir = GetDirectoryArgument(arg)
            else:
                file_path = os.path.abspath(arg)
                directory = os.path.dirname(file_path)
                if not os.path.exists(directory):
                    print("Directory points to the file doesn't exist: %s" % directory)
                    util.Fail()
                if not os.path.isfile(file_path):
                    print("Specified file doesn't exist: %s" % file_path)
                    util.Fail()

        elif arg == "-kp" or arg == "--keep-prev-results":
            options.keepprev = True

        elif arg == "-blfix" or arg == "--fix-branch-length":
            options.fix_branch_length = True

        elif arg == "-si" or arg == "--system-info":
            options.system_info = True
        
        elif arg == "--override":
            options.override = False

        elif arg == "-rs" or arg == "--restart":
            options.restart_from = int(args.pop(0))
        
        elif arg == "-uc" or arg == "--update-cycle":
            options.update_cycle = int(args.pop(0))

        elif arg == "--seed":
            options.seed = int(args.pop(0))

        elif arg == "--mcs-seed":
            options.mcs_seed = int(args.pop(0))

        elif arg == "-iv":
            options.iqtree_version = args.pop(0)

        elif arg == "-una" or arg == "--usr-node-aln":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_node_aln = GetFileArgument(args.pop(0))

        elif arg == "-us" or arg == "--usr-state":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_state = GetFileArgument(args.pop(0))

        elif arg == "-ut" or arg == "--usr-tree":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_tree = GetFileArgument(args.pop(0))

        elif arg == "-uit" or arg == "--usr-iqtree":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_iqtree = GetFileArgument(args.pop(0))

        elif arg == "-utd" or arg == "--usr-mcs-treedir":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_mcs_treeDir = GetDirectoryArgument(args.pop(0))

        elif arg == "-uad" or arg == "--usr-mcs-alndir":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_mcs_alnDir = GetDirectoryArgument(args.pop(0))
        
        elif arg == "-upd" or arg == "--usr-mcs-fadir":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            options.usr_mcs_faDir = GetDirectoryArgument(args.pop(0))

        elif arg == "--outgroups":
            arg_outgroup = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            isfile = False
            isdir = False
            arg = args.pop(0)
            outgroups_list: List[str] = []
            outgroups_dict: Dict[str, Union[str, List[str]]] = {}
            try:
                if os.path.isfile(arg):
                    isfile = True
                    outgroup_path = GetFileArgument(arg)
                    with open(outgroup_path) as reader:
                        for line in reader:
                            line = line.replace("\n", "").strip()
                            outgroups_list.append(line)

                elif os.path.isdir(arg):
                    isdir = True
                    outgroup_dir = GetDirectoryArgument(arg)
                    outgroups_dict = {}
                    for file in os.listdir(outgroup_dir):
                        if fnmatch.fnmatch(file, "*.outgroup*"):
                            outgroup_path = os.path.join(outgroup_dir, file)
                            with open(outgroup_path) as reader:
                                outgroup = []
                                for line in reader:
                                    line = line.replace("\n", "").strip()
                                    outgroup.append(line)
                            outgroups_dict[file.split(".", 1)[0]] = outgroup

                else:
                    if "," in arg:
                        outgroups_list = [og.strip() for og in arg.split(",")]
                    elif "|" in arg:
                        outgroups_list = [og.strip() for og in arg.split("|")]
                    else:
                        outgroups_list = [arg]
            except:
                print("Invalid argument for option %s: %s" % (arg_outgroup, arg))
                print('Valid options are for instance "Nuphar japonica,Nymphaea jamesoniana"\n')
                util.Fail()

            if len(outgroups_list) == 0 and len(outgroups_dict) == 0 and not isfile and not isdir:
                print("Missing option for command line argument from file %s\n" % arg)
                util.Fail()

            if outgroups_dict:      
                options.outgroups = outgroups_dict
            elif outgroups_list:
                options.outgroups = outgroups_list

        elif arg == "-rd":
            if len(args) == 0:
                print("Missing option for command line argument %s" % arg)
                util.Fail()
            options.recDir = GetDirectoryArgument(args.pop(0))

        elif arg == "-st":
            arg_st = arg
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()

            arg = args.pop(0).upper()
            try:
                if len(arg) == len("CODON"):
                    continue
                elif "CODON" in arg and len(arg) > len("CODON"):
                    arg1, arg2 = arg[:5], int(arg[5:])
                    if arg1 == "CODON" and arg2 in [*range(1, 26)]:
                        options.sequence_type = arg
                elif "AA" in arg:
                    options.sequence_type = arg
                
            except:
                print("Invalid argument for option %s: %s" % (arg_st, arg))
                print("Valid options are for instance 'CODON[1-11]' or 'AA'\n")
                print("For more information please refer to http://www.iqtree.org/doc/Substitution-Models#codon-models")
                util.Fail()

        elif arg == "-mbs" or arg == "--mcs-batch-size":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()            

            options.mcs_batch_size = int(args.pop(0))

        elif arg == "-bb":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()            

            options.bootstrap = int(args.pop(0))

        elif arg == "-alrt":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            
            options.sh_alrt = int(args.pop(0))

        elif arg == "-m" or arg == "--evolution-model":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            try:
                options.evolution_model = args.pop(0).upper()
            except:
                print("The substitution models follow the convetion in IQ-TREE2")
                print("For more information please refer to: http://www.iqtree.org/doc/Substitution-Models")
                util.Fail()

        elif arg == "-te":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            
            isfile = False
            isdir = False
            arg = args.pop(0)
            gene_trees = {}
            try:
                if os.path.isfile(arg):
                    isfile = True
                    gene_tree_path = GetFileArgument(arg)
                    file = os.path.basename(gene_tree_path)
                    gene_trees[file.split(".", 1)[0]] = gene_tree_path

                elif os.path.isdir(arg):
                    isdir = True
                    gene_tree_dir = GetDirectoryArgument(arg)
                    for file in os.listdir(gene_tree_dir):
                        if fnmatch.fnmatch(file, "*.tree*") and not fnmatch.fnmatch(file, "*.outgroup*"):
                            gene_tree_path = os.path.join(gene_tree_dir, file)
                            if validate_newick_tree(gene_tree_path):
                                gene_trees[file.split(".", 1)[0]] = gene_tree_path
                            else:
                                gene_trees[file.split(".", 1)[0]] = ""
            except:
                util.Fail()

            if not isfile and not isdir:
                print("Missing option for command line argument from file %s\n" % arg)
                util.Fail()
            options.gene_tree = gene_trees

        elif arg == "--alisim":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.output_prefix = args.pop(0)
        
        elif arg == "-m":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.evolution_model = args.pop(0)  

        elif arg == "--root-seq":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.root_node = args.pop(0)

        elif arg == "--num-alignments":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()   
            options.nalign = int(args.pop(0))     

        elif arg == "-t" or arg == "--nthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)

            try:
                options.nthreads = int(arg)
            except:
                print("Incorrect argument for number of RECUR threads: %s\n" % arg)
                util.Fail()

        elif arg == "-nt" or arg == "--iqtree-nthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)

            try:
                options.iqtree_nthreads = int(arg)
                usr_iqtree_nthread = True
            except:
                print("Incorrect argument for number of IQTREE threads: %s\n" % arg)
                util.Fail()

        elif arg == "-rt" or arg == "--recur-nthreads":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            arg = args.pop(0)
            options.recur_nthreads = int(arg)
            usr_recur_nthread = True

        elif arg == "-fd":
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            
            options.fd_limit = int(args.pop(0)) 

        elif arg == "-n" or arg == "--name":  
            if options.name:
                print("Repeated argument: -n/--name")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail() 
            options.name = args.pop(0)
            while options.name.endswith("/"): options.name = options.name[:-1]
            if any([symbol in options.name for symbol in [" ", "/"]]): 
                print("Invalid symbol for command line argument %s\n" % arg)
                util.Fail()

        elif arg == "-o" or arg == "--output":  
            if resultsDir_nonDefault != None:
                print("Repeated argument: -o/--output")
                util.Fail()
            if len(args) == 0:
                print("Missing option for command line argument %s\n" % arg)
                util.Fail()
            resultsDir_nonDefault = args.pop(0)

            while resultsDir_nonDefault.endswith("/"): 
                resultsDir_nonDefault = resultsDir_nonDefault[:-1]

            resultsDir_nonDefault += "/"
            # if os.path.exists(resultsDir_nonDefault):
            #     print("ERROR: non-default output directory already exists: %s\n" % resultsDir_nonDefault)
            #     util.Fail()

            if " " in resultsDir_nonDefault:
                print("ERROR: non-default output directory cannot include spaces: %s\n" % resultsDir_nonDefault)
                util.Fail()
            checkDirName = resultsDir_nonDefault

            while checkDirName.endswith("/"):
                checkDirName = checkDirName[:-1]

            path, newDir = os.path.split(checkDirName)
            if path != "" and not os.path.exists(path):
                print("ERROR: location '%s' for results directory '%s' does not exist.\n" % (path, newDir))
                util.Fail()

        elif arg == "-h" or arg == "--help":
            helpinfo.PrintHelp()

        elif arg == "-efn" or arg == "--extended-filename":
            options.extended_filename = True

        else:
            print("Unrecognised argument: %s\n" % arg)
            util.Fail()

    if options.nthreads < 1:
        print("ERROR: Number of '-t' RECUR threads cannot be fewer than 1, got %d" % options.nthreads)
        util.Fail()  

    if options.nthreads > mp.cpu_count() and options.iqtree_nthreads == 1:
        print("\nWARNING: Number of threads exceeds the threadshold accepted by RECUR.")
        print(f"{mp.cpu_count()} threads is used instead of {options.nthreads}")
        options.nthreads = mp.cpu_count()

    if not usr_recur_nthread and usr_iqtree_nthread:
        if options.iqtree_nthreads < mp.cpu_count():
            options.recur_nthreads = mp.cpu_count() // options.iqtree_nthreads
        else:
            print("\nWARNNING: IQ-TREE2 threads exceeds the threadshold accepted by RECUR.")
            print("Adjusting the threads as follows: ")
            options.recur_nthreads, options.iqtree_nthreads = 1, mp.cpu_count()
            print(f"IQ-TREE2 threads: {options.iqtree_nthreads}")

    if usr_iqtree_nthread and usr_recur_nthread:
        if options.recur_nthreads * options.iqtree_nthreads > mp.cpu_count():
            print("\nWARNING: The combination of the number of RECUR threads and IQ-TREE2 threads exceeds the threadshold accepted by RECUR.")
        
        if options.iqtree_nthreads <= mp.cpu_count():
            options.recur_nthreads = mp.cpu_count() // options.iqtree_nthreads
        
        elif options.iqtree_nthreads > mp.cpu_count():
            options.recur_nthreads, options.iqtree_nthreads = 1, mp.cpu_count()
        
        print("Adjusting the threads as follows: ")
        print(f"RECUR threads: {options.recur_nthreads}, IQ-TREE2 threads: {options.iqtree_nthreads}")
        print("\nRECOMMENDATION: To ensure reproducibility and optimal performance, configure the number of threads used by IQ-TREE2 to be either 1 or to match the number of threads used by RECUR.")
     
    if not options.qStartFromMSA:
        print("ERROR: Please specify the input directory for RECUR using the option: '-f'.")
        util.Fail()

    if resultsDir_nonDefault:
        options.recDir = resultsDir_nonDefault

    return options, alnDir, alnPath, resultsDir_nonDefault

