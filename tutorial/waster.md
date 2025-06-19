# Without-Alignment/Assembly Species Tree EstimatoR † (WASTER)

[<img src="../misc/WASTER.png" width="500"/>](../misc/WASTER.png)

WASTER is a coalesence-aware ***de novo*** species tree inference tool, which means it can take as inputs raw reads in FASTQ format.
Paticularly, WASTER can accurately infer species tree even from Illumina reads with only ***1.5X depth***.
WASTER infers the species tree by first calling SNPs from reads/assembies and then invoking CASTER to reconstruct the species tree from the SNPs.

## Publication

Chao Zhang, Rasmus Nielsen, WASTER: Practical de novo phylogenomics from low-coverage short reads, bioRxiv (2025) https://doi.org/10.1101/2025.01.20.633983


# Announcements

## Integrated in Phylosuite (NEW)

Many ASTER tools have been integrated in [PhyloSuite](http://phylosuite.jushengwu.com/), an integrated and scalable desktop platform for streamlined molecular sequence data management and evolutionary phylogenetics studies.

## GUI for Windows users

Please check out our software with GUI. Simply download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip), extract the contents, enter `exe` folder, and click `aster-gui.exe`. 

## Bug Reports

Contact ``chao.zhang@sund.ku.dk``, [``aster-users@googlegroups.com``](https://groups.google.com/forum/#!forum/aster-users), or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this TUTORIAL file
- Forums (feel free to ask questions or ask for help running ASTER):
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER issues page](https://github.com/chaoszhang/ASTER/issues)
  - QQ group: 130635706

# INSTALLATION
For most users, installing ASTER is ***very*** easy!
Download using one of two approaches:
  - You simply need to download the zip file for [Windows](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip)/[MacOS](https://github.com/chaoszhang/ASTER/archive/refs/heads/MacOS.zip)/[Linux](https://github.com/chaoszhang/ASTER/archive/refs/heads/Linux.zip) and extract the contents to a folder of your choice.
  - Alternatively, you can clone the [github repository](https://github.com/chaoszhang/ASTER.git) and checkout the branch named Windows/MacOS/Linux.

Binary files should be in the `exe` folder for Windows or `bin` folder otherwise. If you are lucky, these may just work as is and you may not need to build at all.

## For Linux/Unix/WSL users
1. In terminal, `cd` into the downloaded directory and run `make`.
  - If you see `*** Installation complete! ***` then you are done!
  - If you see `Command 'g++' not found` then before rerunning `make`,
    - Debian (Ubuntu) users try
      ```
      sudo apt update
      sudo apt install g++
      ```
    - CentOS (RedHat) users try
      ```
      sudo yum update
      sudo yum install gcc-c++
      ```
    - Unix (MacOS) users should be prompted for installing `g++` and please click "install". If no prompt, try `g++`.
  - If you see "error" when running `make`, please try `make waster` instead and file a bug report.
2. Binary files should be in the `bin` folder.

## For Windows users
- [Executables](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip) for x86-64 are available in `exe` folder and it is **very likely** that they already work.
- [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) is HIGHLY recommanded if you need to install on your own! Please follow instructions in "For Linux/Unix/WSL users" section.
- To compile windows excutables:
  1. Download [MinGW](https://sourceforge.net/projects/mingw-w64/) and install ***posix*** version for your architecture (eg. x86-64)
  2. Add path to `bin` folder of MinGW to [system environment variable `PATH`](https://www.google.com/search?q=Edit+the+system+environment+variables+windows)
  3. Double click `make.bat` inside the downloaded directory

### GUI for Windows users (NEW)

Please check out our software with GUI. Simply download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip), extract the contents, enter `exe` folder, and click `aster-gui.exe`. 

# STOP!
If you have ***overlapping*** paired-end sequencing reads, please make sure you have merged them using [`BBMerge`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) or alike.
This will improve the accuracy of WASTER on biological datasets!
If you have non-overlapping paired-end sequencing reads, great!
With the same sequencing depth, WASTER prefers non-overlapping reads.

A detailed [walkthrough](../misc/waster-linux-walkthrough.md) is also available (English/中文).

# INPUT
The input is a text file containing a list of species names each followed by a Fasta/Fastq file (one species name and one file per line). For example:
```
speciesA	example/waster/filename1.fa
speciesB	example/waster/filename2.fa
speciesC	example/waster/filename3.fq
speciesD	example/waster/filename4.fq
```
Make sure input files are ***not zipped***.


When multiple individuals from the same species are available, you need to let WASTER to know that they are from the same species.
  1. You can give multiple individuals from the same species the same name in the input list.
```
speciesA	example/waster/filename1.fa
speciesA	example/waster/filename2.fa
speciesB	example/waster/filename3.fq
speciesB	example/waster/filename4.fq
...
```
***WARNING: Never use this method to input pair-end reads files!!!***
If your pair-end reads overlap, then merge them using [`BBMerge`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) or alike; otherwise `cat` them into one file.

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch supports measured as local bootstrap support (>95.0 means good)
* It can also annotate branches with other quantities, such as quartet scores and local bootstraps for all three topologies.

# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER and issue the following command:

```
bin/waster
```

This will give you a list of options available. If you are using Windows, please replace `bin/waster` with `.\exe\waster.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
bin/waster INPUT_FILE
```
or
```
bin/waster -i INPUT_FILE
```

In the first case, INPUT_FILE is ***hard-coded*** to be the ***last argument*** for backward compatibility. 

For example if you want to run `waster` with input `example/waster/input_list.txt`, then run

```
bin/waster example/waster/input_list.txt
```
or
```
bin/waster -i example/waster/input_list.txt
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
bin/waster -o OUTPUT_FILE INPUT_FILE
```
or
```
bin/waster -i INPUT_FILE -o OUTPUT_FILE
```

With `-i INPUT_FILE` option, the order does not matter anymore. For brevity, from here on we will not demonstrate `-i INPUT_FILE` cases.

To save the logs (**also recommended**), run:

```
bin/waster -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

For example, you can run

```
bin/waster -o example/waster/input_list.txt.stree example/waster/input_list.txt 2>example/waster/input_list.txt.log
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
bin/waster -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

ASTER also allows rooting at an given outgroup:

```
bin/waster --root YOUR_OUTGROUP INPUT_FILE
```

***Notice:*** By default, WASTER will create a 32 GB look-up table to find SNPs. So, if you are running on a machine with <64 GB memory, you need to shrink the look-up table size using `-k 8` or `-k 7` even for the example run.
Using `-k 8` requires about 4 GB memory and `-k 7` only requires <1 GB memory.
Try the following example run:
```
bin/waster-site -k 7 example/waster/input_list.txt
```

## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
bin/waster -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to run with more rounds of placement for ensured optimality, then you can run with
```
bin/waster -r 16 -s 16 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```
or simply
```
bin/waster -R -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/waster -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`. If want to score a species tree or you want to place only ***one*** taxon onto the tree, you can use

```
bin/waster -r 1 -s 0 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```
or simply,
```
bin/waster -C -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/waster -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```
