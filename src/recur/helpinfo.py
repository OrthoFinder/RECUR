
from recur.citation import print_citation
from recur import recur_iqtree_nthreads, iqtree_nthreads, recur_nthreads

def PrintHelp():  

    print("")
    print("SIMPLE USAGE:") 
    print("Run full RECUR analysis on a protein or codon alignment <dir/file>")
    print("recur [options] -f <dir/file> --outgroups <outgroup_species/dir/file> -st <AA|CODON>")
    
    print("")
    print("OPTIONS:")
    print("-f <dir/file>                Protein or codon alignment in FASTA format [Required]")
    print("-s <str>                     <AA|CODON> [Required][Default: CODON1]")
    print("--outgroups <dir/file/str>   List of outgroup species [Required]")
    print("--num-alignments <int>       Number of simulated alignments for p-value estimation [Default: 1000]")
    print("-te <dir/file>               Complete constraint tree [Default: Estimated from alignment]")
    print("-m <str>                     Model of sequence evolution [Default: estimated from alignment]")
    print(f"-nt <int>                    Number of threads provided to IQ-TREE2 [Default: 1 (without alrt); {iqtree_nthreads} (with alrt)]")
    print(f"-rt <int>                    Number of threads used for RECUR run on IQ-TREE2 [Default: {recur_iqtree_nthreads}]")
    print(f"-t <int>                     Number of threads used for RECUR internal processing [Default: {recur_nthreads}]")
    print("--seed <int>                 Random starting see number [Default: 8]")
    print("-o <txt>                     Results directory [Default: same directory as MSA files]")
    print("-pbs <int>                   Batch size for substituion analysis on protein alignments of real phylogeny [Default: None]") 
    print("-mbs <int>                   Batch size for substituion analysis on the Monte Carlo Simulated alignments [Default: None]")
    print("-iv <str>                    IQ-TREE2 path [Default: local]")
    
    
    print("")
    print("LICENSE:")
    print(" Distributed under the GNU General Public License (GPLv3). See License.md")
    print(print_citation) 
