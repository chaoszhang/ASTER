# Weighted ASTRAL (wASTRAL)
Weighted ASTRAL program has the following modes for different weighting schemes:
1. Weighted ASTRAL by Branch Support (mode 2)
2. Weighted ASTRAL by Branch Length (mode 3)
3. Weighted ASTRAL - Hybrid (default)

Weighted ASTRAL series introduce threshold-free weighting schemes for the quartet-based species tree inference, the metric used in the popular method ASTRAL. By reducing the impact of quartets with low support or long terminal branches (or both), weighting provides stronger theoretical guarantees and better empirical performance than the unweighted ASTRAL. Our results show that weighted ASTRAL improves the utility of summary methods and can reduce the incongruence often observed across analytical pipelines.

## Publication

[1] Chao Zhang, Siavash Mirarab, Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees, Molecular Biology and Evolution, 2022, msac215, https://doi.org/10.1093/molbev/msac215

(OPTIONAL) [2] Chao Zhang, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).

### Example of usage

We obtained the species tree from gene trees using wASTRAL v1.22.3.7 [1]. 
(OPTIONAL) This method reimplements ASTRAL-III [2] and takes into account phylogenetic uncertainty by intergrating signals from branch length and branch support in gene trees.


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
  - If you see "error" when running `make`, please try `make wastral` instead and file a bug report.
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

# INPUT
* The input gene trees are in the Newick format
* The input trees can have missing taxa, polytomies (unresolved branches), and multiple individuals/genes per species.
* When individuals/genes from the same species are available, you can ask ASTRAL to force them to be together in the species tree. You can do this in two ways.
  1. You can give multiple individuals/genes from the same species the same name in the input gene trees.
  2. OR, a mapping file needs to be provided using the `-a` option.
```
individual_A1 species_name_A
individual_A2 species_name_A
individual_B1 species_name_B
individual_B2 species_name_B
individual_B3 species_name_B
...
```
  Or
```
gene_A1 species_name_A
gene_A2 species_name_A
gene_B1 species_name_B
gene_B2 species_name_B
gene_B3 species_name_B
...
```

1. **Weighted ASTRAL by Branch Support (mode 2)**: Non-root interal node labels must be a non-negative number being support. Eg. `((A,B)100,(C,D)0);` or `((A:1,B:1)1.0:1,(C:1,D:1)0.333:0);`.
2. **Weighted ASTRAL by Branch Length (mode 3)**: Non-root labels must have branch lengths after `:`. Eg. `((A,B):1,(C,D):0);` or `((A:1,B:1)1.0:1,(C:1,D:1)0.333:0);`.
3. **Weighted ASTRAL - Hybrid (default)**: Non-root interal node labels must be a non-negative number being support before `:` and non-root labels must have branch lengths after `:`. Eg. `((A:1,B:1)100:1,(C:1,D:1)0:0);` or `((A:1,B:1)1.0:1,(C:1,D:1)0.333:0);`.

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch lengths in coalescent units for astral-weighted and in **combined (coalensecent + 2 * substitution) units** (only for internal branches)
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.

The weighted ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).

# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER and issue the following command:

```
bin/wastral
```

This will give you a list of options available. If you are using Windows, please replace `bin/wastral` with `.\exe\wastral.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
bin/wastral INPUT_FILE
```
or
```
bin/wastral -i INPUT_FILE
```

In the first case, INPUT_FILE is ***hard-coded*** to be the ***last argument*** for backward compatibility. 

For example if you want to run `wastral` with input `example/genetree.nw`, then run

```
bin/wastral example/genetree.nw
```
or
```
bin/wastral -i example/genetree.nw
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
bin/wastral -o OUTPUT_FILE INPUT_FILE
```
or
```
bin/wastral -i INPUT_FILE -o OUTPUT_FILE
```

With `-i INPUT_FILE` option, the order does not matter anymore. For brevity, from here on we will not demonstrate `-i INPUT_FILE` cases.

To save the logs (**also recommended**), run:

```
bin/wastral -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

For example, you can run

```
bin/wastral -o example/genetree.nw.stree example/genetree.nw 2>example/genetree.nw.log
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
bin/wastral -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

ASTER also allows rooting at an given outgroup:

```
bin/wastral --root YOUR_OUTGROUP INPUT_FILE
```

By default, wASTRAL assumes multiple individuals/alleles from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/wastral -a example/genetree.map example/genetree.nw
```

When your dataset has no more than 50 species and no more than 500 genes, you may want to run with more rounds using `-R` (see below). 

## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
bin/wastral -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to run with more rounds of placement for ensured optimality, then you can run with
```
bin/wastral -r 16 -s 16 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```
or simply
```
bin/wastral -R -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/wastral -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`. If want to score a species tree or you want to place only ***one*** taxon onto the tree, you can use

```
bin/wastral -r 1 -s 0 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```
or simply,
```
bin/wastral -C -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/wastral -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```

Add `-u 0` before `INPUT_FILE` if you want to compute species tree topology only; Add `-u 2` before `INPUT_FILE` if you support and local-PP for all three resolutions of each branch.

```
bin/wastral -u 0 -o OUTPUT_FILE INPUT_FILE
bin/wastral -u 2 -o OUTPUT_FILE INPUT_FILE
```

Species tree with more than **2000** taxa may cause **floating point underflow or precision issue**. Use the following command instead:

```
make wastral_precise
bin/wastral_precise -o OUTPUT_FILE INPUT_FILE
```

***Notice:*** For hybrid weighting (default) and weighting by support (mode 2), wASTRAL by default will automatically detect support type.
You may also specify max (`-x`) and min (`-n`) of support value to the program.

```
bin/wastral -x MAX_SUPPORT -n MIN_SUPPORT INPUT_FILE
```

For **Bootstrap** support, any of the following commands works:

```
bin/wastral -S INPUT_FILE
bin/wastral -x 100 -n 0 INPUT_FILE
```

For **local Baysian** support, `-x 1 -n 0.333` is recommended or `-B` for short:

```
bin/wastral -B INPUT_FILE
```
or
```
bin/wastral -x 1 -n 0.333 INPUT_FILE
```

For other **probability & likelihood** support, `-x 1 -n 0` may be more reasonable or `-L` for short:

```
bin/wastral -L INPUT_FILE
```
or
```
bin/wastral -x 1 -n 0 INPUT_FILE
```

By default, wASTRAL assumes multiple individuals/alleles from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/wastral -a example/genetree.map example/genetree.nw
```

In case you want to run Weighted ASTRAL by Branch Support, you can do:

```
bin/wastral --mode 2 -o OUTPUT_FILE INPUT_FILE
```

In case you want to run Weighted ASTRAL by Branch Length, you can do:

```
bin/wastral --mode 3 -o OUTPUT_FILE INPUT_FILE
```

