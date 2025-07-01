
from recur.utils.process_args import iqtree_nthreads, recur_nthreads
from recur.citation import print_citation
try:
    from rich import print
except ImportError:
    ...
from rich.table import Table
from rich.console import Console


width = 28

def PrintHelp(other_options: bool = False) -> None:

    console = Console()
    table_options = Table(show_header=False, box=None, expand=False)
    table_options.add_column("Option", justify="left", width=width, no_wrap=False, overflow="fold")
    table_options.add_column("Description", justify="left", overflow="fold")


    other_options_table = Table(show_header=False, box=None, expand=False)
    other_options_table.add_column("Option", justify="left", width=width, no_wrap=False, overflow="fold")
    other_options_table.add_column("Description", justify="left", overflow="fold")

    print("")
    print("SIMPLE USAGE:")
    print("Run full [dark_goldenrod]RECUR[/dark_goldenrod] analysis on a protein or codon alignment <[bright_magenta]dir/file[/bright_magenta]>")
    print("recur [options] -f <[bright_magenta]dir/file[/bright_magenta]> --outgroups <[bright_magenta]dir/file/str[/bright_magenta]> -st <[dark_cyan]AA|CODON[/dark_cyan]>")

    table_options.add_row(
        "-f <[bright_magenta]dir/file[/bright_magenta]>",
        f"Protein or codon alignment in [red1]FASTA[/red1] format [orange3][Required][/orange3]"
    )

    table_options.add_row(
        "-st <[bright_magenta]str[/bright_magenta]>",
        f"<[dark_cyan]AA|CODON[/dark_cyan]> [orange3][Required][/orange3][Default: [dark_cyan]AA[/dark_cyan]]"
    )

    table_options.add_row(
        "--outgroups <[bright_magenta]dir/file/str[/bright_magenta]>",
        f"List of outgroup sequences [orange3][Required][/orange3]"
    )

    table_options.add_row(
        "--num-alignments <[bright_magenta]int[/bright_magenta]>",
        f"Number of simulated alignments for p-value estimation [Default: [deep_sky_blue2]1000[/deep_sky_blue2]]"
    )

    table_options.add_row(
        "-te <[bright_magenta]dir/file[/bright_magenta]>",
        f"Complete constraint tree [Default: [dark_cyan]estimated from alignment[/dark_cyan]]"
    )

    table_options.add_row(
        "-m <[bright_magenta]str[/bright_magenta]>",
        f"Model of sequence evolution [Default: [dark_cyan]estimated from alignment[/dark_cyan]]"
    )


    table_options.add_row(
        "-nt <[bright_magenta]int[/bright_magenta]>",
        f"Number of threads provided to IQ-TREE [Default: [deep_sky_blue2]1[/deep_sky_blue2] (without alrt); [deep_sky_blue2]{iqtree_nthreads}[/deep_sky_blue2] (with alrt)]"
    )



    table_options.add_row(
        "--seed <[bright_magenta]int[/bright_magenta]>",
        f"Random starting see number [Default: [deep_sky_blue2]8[/deep_sky_blue2]]"
    )

    table_options.add_row(
        "-o <[bright_magenta]txt[/bright_magenta]>",
        f"Results directory [Default: [dark_cyan]same directory as MSA files[/dark_cyan]]"
    )
    table_options.add_row(
        "-uc <[bright_magenta]int[/bright_magenta]>",
        f"Update cycle used in progress bar [Default: [dark_cyan]no progress bar[/dark_cyan]]"
    )

    table_options.add_row(
        "-bs <[bright_magenta]int[/bright_magenta]>",
        f"Batch size used in subsitution analysis of the Monte Carlo Simulated sequences [Default: [dark_cyan]no batch processing[/dark_cyan]]"
    )

    table_options.add_row(
        "-iv <[bright_magenta]str[/bright_magenta]>",
        f"IQ-TREE version. [Default: [dark_cyan]iqtree2[/dark_cyan]]"
    )

    table_options.add_row(
        "-sl <[bright_magenta]float[/bright_magenta]>",
        f"Significance level. [Default: [deep_sky_blue2]0.05[/deep_sky_blue2]]"
    )

    table_options.add_row(
        "-pam <[bright_magenta]str[/bright_magenta]>",
        f"P-Value adjustment method. Available methods: bonferroni, holm, fdr_bh, fdr_by, fdr_tsbh. [Default: [dark_cyan]None[/dark_cyan]]"
    )

    table_options.add_row(
        "-blfix",
        f"Fix branch lengths of tree. [Default: [dark_cyan]False[/dark_cyan]]"
    )
    table_options.add_row(
        "-nbt",
        f"Branch test control. [Default: [dark_cyan]True[/dark_cyan]]"
    )

    table_options.add_row(
        "-ps",
        f"Return the P-values statistics in the recurrence list. [Default: [dark_cyan]False[/dark_cyan]]"
    )

    print("")
    console.print("[bold]OPTIONS:[/bold]")
    console.print(table_options)
    console.print()



    if other_options:
        
        other_options_table.add_row(
            "-rs <[bright_magenta]int[/bright_magenta]>",
            f"Restart [dark_goldenrod]RECUR[/dark_goldenrod] is a specified step. [Default: [deep_sky_blue2]1[/deep_sky_blue2]]"
        )

        other_options_table.add_row(
            "-rl <[bright_magenta]int[/bright_magenta]>",
            f"The number MCS that [dark_goldenrod]RECUR[/dark_goldenrod] can handle in a single stage. [Default: [dark_cyan]1000[/dark_cyan]]"
        )

        other_options_table.add_row(
            "-ret <[bright_magenta]float[/bright_magenta]>",
            f"Relative for the P-values which affects the number of MCS. [Default: [deep_sky_blue2]0.1[/deep_sky_blue2]]"
        )

        other_options_table.add_row(
            "-mce",
            f"Monte Carlo Error control. [Default: [dark_cyan]False[/dark_cyan]]"
        )
        other_options_table.add_row(
            "-gc",
            f"Grid cushion which will double the computed size of MCS. [Default: [dark_cyan]False[/dark_cyan]]"
        )

        other_options_table.add_row(
            "-ps",
            f"Return the P-values statistics. [Default: [dark_cyan]False[/dark_cyan]]"
        )

        other_options_table.add_row(
            "-jr",
            f"Return only the reccurence list for the real phylogeny. [Default: [dark_cyan]False[/dark_cyan]]"
        )
        other_options_table.add_row(
            "-ds",
            f"Removes the MC simulated MSA files after extracting the recurrence count. [Default: [dark_cyan]False[/dark_cyan]]"
        )
        other_options_table.add_row(
            "-ms",
            f"Run [dark_goldenrod]RECUR[/dark_goldenrod] in multi stages. [Default: [dark_cyan]False[/dark_cyan]]"
        )

        other_options_table.add_row(
            "-cr",
            f"Compute recurrence based on the MC simulated MSA recurrence count file. Used in multi-stage [dark_goldenrod]RECUR[/dark_goldenrod]. [Default: [dark_cyan]False[/dark_cyan]]"
        )

        print("")
        console.print("[bold]OTHER OPTIONS:[/bold]")
        console.print(other_options_table)
        console.print()

    # print("-f <dir/file>                Protein or codon alignment in FASTA format [Required]")
    # print("-st <str>                     <AA|CODON> [Required][Default: CODON1]")
    # print("--outgroups <dir/file/str>   List of outgroup sequences [Required]")
    # print("--num-alignments <int>       Number of simulated alignments for p-value estimation [Default: 1000]")
    # print("-te <dir/file>               Complete constraint tree [Default: estimated from alignment]")
    # print("-m <str>                     Model of sequence evolution [Default: estimated from alignment]")
    # print("-blfix                       Fix branch lengths of tree. [Default: False]")
    # print(f"-nt <int>                    Number of threads provided to IQ-TREE [Default: 1 (without alrt); {iqtree_nthreads} (with alrt)]")
    # print(f"-t <int>                     Number of threads used for RECUR internal processing [Default: {recur_nthreads}]")
    # print("--seed <int>                 Random starting see number [Default: 8]")
    # print("-o <txt>                     Results directory [Default: same directory as MSA files]")
    # print("-uc <int>                    Update cycle used in progress bar [Default: no progress bar]")
    # print("-bs <int>                    Batch size used in subsitution analysis of the Monte Carlo Simulated sequences [Default: no batch processing]")


    print("")
    print("[bold]LICENSE:[/bold]")
    print(" Distributed under the [dodger_blue1]GNU General Public License (GPLv3)[/dodger_blue1]. See License.md")
    print(print_citation)
    table_options.add_row(
        "-t <[bright_magenta]int[/bright_magenta]>",
        f"Number of threads used for RECUR internal processing [Default: [deep_sky_blue2]{recur_nthreads}[/deep_sky_blue2]]"
    )