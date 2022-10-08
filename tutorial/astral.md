# Accurate Species Tree ALgorithm (ASTRAL)
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of tripartitions in the species tree comes from a predefined set of tripartitions.
ASTER re-implements [ASTRAL](https://github.com/smirarab/ASTRAL) as a complementary to the original ASTRAL on datasets for which the original ASTRAL is not suitable (e.g. large datasets, multi-individual, and super-tree construction).

Warning: ASTER implementation may be **slower** and even **less accurate** than ASTRAL-III when input gene trees has fewer than 50 species and 500 genes. Please choose wisely!
As a supplementary to ASTRAL-III, ASTER lacks of many features of ASTRAL-III (e.g. computing support). You can work around by first computing optimal tree with ASTER and use the ASTER output tree as `-q` option to ASTRAL-III for annotation. 

## Publication

Chao Zhang, Siavash Mirarab, Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees, Molecular Biology and Evolution, 2022, msac215, https://doi.org/10.1093/molbev/msac215

Zhang, Chao, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).

## Bug Reports

Contact ``aster-users@googlegroups.com`` or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this TUTORIAL file
- [README/astral-pro.md](README/astral-pro.md) for ASTRAL and ASTRAL-Pro; [README/wastral.md](README/wastral.md) for weighted ASTRAL series; [README/asterisk.md](README/asterisk.md) for ASTERISK series
- Forums:
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER issues page](https://github.com/chaoszhang/ASTER/issues)

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
    - Unix (MacOS) users should be prompted for installing `g++` and please click "install". If no prompt, try `g++`
2. Binary files should be in the `bin` folder.

## For Windows users
- Executables for x86-64 are available in `exe` folder and it is **very likely** that they already work.
- [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) is HIGHLY recommanded if you need to install on your own! Please follow instructions in "For Linux/Unix/WSL users" section.
- To compile windows excutables:
  1. Download [MinGW](https://sourceforge.net/projects/mingw-w64/) and install ***posix*** version for your architecture (eg. x86-64)
  2. Add path to `bin` folder of MinGW to [system environment variable `PATH`](https://www.google.com/search?q=Edit+the+system+environment+variables+windows)
  3. Double click `make.bat` inside the downloaded directory

### GUI for Windows users (NEW)

Please check out our software with GUI. Simply download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip), extract the contents, enter `exe` folder, and click `aster-gui.exe`. 

# INPUT
* The input trees can have missing taxa, polytomies (unresolved branches), and multiple individuals per species.
* When individuals genes from the same species are available, you can ask ASTRAL to force them to be together in the species tree. You can do this in two ways.
  1. You can give multiple individuals from the same species the same name in the input gene trees (e.g., `((species_name_A,species_name_B),(species_name_A,species_name_C));`).
  2. OR, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per genes, and each line needs to be in the following formats (e.g., for gene trees like `((individual_A1,individual_B1),(individual_A2,individual_C1));`):
```
individual_A1 species_name_A
individual_A2 species_name_A
individual_B1 species_name_B
individual_B2 species_name_B
individual_B3 species_name_B
...
```

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch lengths in coalescent units (only for internal branches)
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.

The ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).

# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER and issue the following command:

```
bin/astral
```

This will give you a list of options available. If you are using Windows, please replace `bin/astral` with `.\exe\astral.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
bin/astral INPUT_FILE
```
or
```
bin/astral -i INPUT_FILE
```

In the first case, INPUT_FILE is ***hard-coded*** to be the ***last argument*** for backward compatibility. 

For example if you want to run `astral` with input `example/genetree.nw`, then run

```
bin/astral example/genetree.nw
```
or
```
bin/astral -i example/genetree.nw
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
bin/astral -o OUTPUT_FILE INPUT_FILE
```
or
```
bin/astral -i INPUT_FILE -o OUTPUT_FILE
```

With `-i INPUT_FILE` option, the order does not matter anymore. For brevity, from here on we will not demonstrate `-i INPUT_FILE` cases.

To save the logs (**also recommended**), run:

```
bin/astral -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

For example, you can run

```
bin/astral -o example/genetree.nw.stree example/genetree.nw 2>example/genetree.nw.log
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
bin/astral -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

By default, ASTRAL assumes multiple individuals from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/astral -a example/genetree.map example/genetree.nw
```

When your dataset has no more than 50 species and no more than 500 genes, you may want to run with more rounds using `-R` (see below). 

## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
bin/astral -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to run with more rounds of placement for ensured optimality, then you can run with
```
bin/astral -r 16 -s 16 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```
or simply
```
bin/astral -R -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/astral -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`. If want to score a species tree or you want to place only ***one*** taxon onto the tree, you can use

```
bin/astral -r 1 -s 0 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```
or simply,
```
bin/astral -C -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/astral -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```

Species tree with more than **5000** taxa may cause **overflow**. Use the following command instead:

```
bin/astral_int128 -o OUTPUT_FILE INPUT_FILE
```
