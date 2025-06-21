import os
from typing import Dict, List, Optional

from recur.utils import parallel_task_manager, util

# asr_command_str = "iqtree2 -s alignment_file -redo -T iqtree_nthreads -m evolution_model -pre path_to_output --seed phy_seed -safe"
# Monte_Carlo_simulation_command_str = "iqtree2 --alisim output_prefix -T iqtree_nthreads -m best_evolution_model -te gene_tree --keep-seq-order --root-seq root_node --num-alignments nalign --seed mcs_seed --write-all --out-format fasta -safe"

asr_command = lambda iqtree_version: f"{iqtree_version} -s alignment_file -redo -T iqtree_nthreads -m evolution_model -pre path_to_output --seed phy_seed -safe"
Monte_Carlo_simulation_command = lambda iqtree_version: f"{iqtree_version} --alisim output_prefix -T iqtree_nthreads -m best_evolution_model -te gene_tree --keep-seq-order --root-seq root_node --num-alignments nalign --seed mcs_seed --write-all --out-format fasta -safe"

def GetGeneTreeBuildCommands(
        alignment_file_list: list[str],
        output_dir: str,
        evolution_model: str,
        iqtree_nthreads: int,
        phy_seed: int = 8,
        asr: bool = False,
        output_prefix: Optional[str] = None,
        sequence_type: Optional[str] = None,
        gene_tree: Optional[str] = None,
        bootstrap: Optional[int] = None,
        sh_alrt: Optional[int] = None,
        branch_test: bool = True,
        fix_branch_length: bool = False,
        iqtree_version: str = "iqtree2",
        iqtree_cmd_dict: Dict[str, Dict[str, str]] = {}
    ) -> list[str]:

    commands = []
    for alignment_file in alignment_file_list:
        identifier = os.path.basename(alignment_file).rsplit(".", 1)[0]
        if not output_prefix:
            pre = output_dir + identifier
        else:
            pre = output_prefix + identifier

        command = GetGeneTreeBuildCommand(
            alignment_file,
            iqtree_nthreads,
            evolution_model,
            pre,
            phy_seed,
            asr=asr,
            sequence_type=sequence_type,
            gene_tree=gene_tree,
            bootstrap=bootstrap,
            sh_alrt=sh_alrt,
            branch_test=branch_test,
            fix_branch_length=fix_branch_length,
            iqtree_version=iqtree_version,
            iqtree_cmd_dict=iqtree_cmd_dict,
        )

        commands.append(command)
    return commands

def GetGeneTreeBuildCommand(
        alignment_file: str,
        iqtree_nthreads: int,
        evolution_model: str,
        path_to_output: str,
        phy_seed: int,
        asr: bool = False,
        sequence_type: Optional[str] = None,
        gene_tree: Optional[str] = None,
        bootstrap: Optional[int] = None,
        sh_alrt: Optional[int] = None,
        branch_test: bool = True,
        fix_branch_length: bool = False,
        iqtree_version: str = "iqtree2",
        iqtree_cmd_dict: Dict[str, Dict[str, str]] = {}
    ) -> str:
    
    if len(iqtree_cmd_dict) == 0:
        asr_command_str = asr_command(iqtree_version)
    else:
        try:
            asr_command_str = iqtree_cmd_dict[iqtree_version]["asr_cmd"]
        except:
            print(f"Unable to run ancestral state reconstruction.")
            print(f"Cannot find {iqtree_version} in the configuration file.")
            util.Fail()

    command = asr_command_str.\
                        replace("alignment_file", alignment_file).\
                        replace("phy_seed", str(phy_seed)).\
                        replace("iqtree_nthreads", str(iqtree_nthreads)).\
                        replace("evolution_model", evolution_model).\
                        replace("path_to_output", path_to_output).split()  # replace("nthreads", str(nthreads)).\

    if sequence_type:
        command.extend(["-st", sequence_type])

    if asr:
        command.append("-asr")

    if gene_tree:
        command.extend(["-te", gene_tree]) # if this is provided, we don't need -bb and -alrt
        bootstrap = None
        sh_alrt = None
        if "-bb" in command:
            command.remove("-bb")
        if "-alrt" in command:
            command.remove("-alrt")

    if not gene_tree and branch_test: 
        if bootstrap:
            command.extend(["-bb", str(bootstrap)])
        if sh_alrt:
            command.extend(["-alrt", str(sh_alrt)])

    if fix_branch_length:
        command.append("-blfix")

    return " ".join(command)

def GetMCsimulationCommand(
        output_prefix: str,
        iqtree_nthreads: int,
        mcs_seed: int,
        best_evolution_model: str,
        gene_tree: str,
        root_node: str,
        nalign: int,
        iqtree_version: str = "iqtree2",
        iqtree_cmd_dict: Dict[str, Dict[str, str]] = {}
    ) -> List[str]:
    
    if len(iqtree_cmd_dict) == 0:
        Monte_Carlo_simulation_command_str = Monte_Carlo_simulation_command(iqtree_version)
    else:
        try:
            Monte_Carlo_simulation_command_str = iqtree_cmd_dict[iqtree_version]["alisim_cmd"]
        except:
            print(f"Unable to run Monte Carlo Simulation.")
            print(f"Cannot find {iqtree_version} in the configuration file.")
            util.Fail()
    
    command = Monte_Carlo_simulation_command_str.\
                                     replace("output_prefix", output_prefix).\
                                     replace("iqtree_nthreads", str(iqtree_nthreads)).\
                                     replace("mcs_seed", str(mcs_seed)).\
                                     replace("best_evolution_model", best_evolution_model).\
                                     replace("gene_tree", gene_tree).\
                                     replace("root_node", root_node).\
                                     replace("nalign", str(nalign))

    return [command]


def RunCommand(
        commands: list[str],
        fileDir: str,
        env: Optional[Dict[str, str]] = None,
        nthreads: int = 1,
        delete_files: bool = False,
        files_to_keep: Optional[List[str]] = None,
        files_to_remove: Optional[List[str]] = None,
        fd_limit: Optional[int] = None,
    ) -> None:

    parallel_task_manager.RunParallelCommands(
        nthreads, 
        commands, 
        fileDir,
        env=env,
        delete_files=delete_files,
        files_to_keep=files_to_keep,
        files_to_remove=files_to_remove,
        q_print_on_error=True,
        q_always_print_stderr=False,
        fd_limit=fd_limit,
    )
