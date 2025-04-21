# Accurate Species Tree ALgorithm for PaRalogs and Orthologs (ASTRAL-Pro3)
ASTRAL-Pro stands for ASTRAL for PaRalogs and Orthologs. ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees and is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL-pro extends ASTRAL to allow multi-copy genes. ASTRAL-pro finds the species tree that has the maximum number of shared induced quartet tree equivalent classes with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. Please see the paper below for the definition of the PL-quartet scores, which is what ASTRAL-Pro optimizes. We refer to the tool both as A-Pro and ASTRAL-Pro. 

ASTRAL-Pro3 re-implements [ASTRAL-Pro](https://github.com/chaoszhang/A-pro) in an equally accurate yet **faster**, and **easier to install** and **lower memory consumption** way.
ASTRAL-Pro3 also integrates [CASTLES-Pro](https://github.com/ytabatabaee/CASTLES) and thus computes terminal and internal branch lengths in substitution-per-site units.

## Publication

[1] Chao Zhang, Siavash Mirarab, ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees, Bioinformatics, 2022, btac620, https://doi.org/10.1093/bioinformatics/btac620

[2] Chao Zhang, Celine Scornavacca, Erin K Molloy, Siavash Mirarab, ASTRAL-Pro: Quartet-Based Species-Tree Inference despite Paralogy, Molecular Biology and Evolution, Volume 37, Issue 11, November 2020, Pages 3292–3307, https://doi.org/10.1093/molbev/msaa139

[3] Yasamin Tabatabaee, Chao Zhang, Tandy Warnow, Siavash Mirarab, Phylogenomic branch length estimation using quartets, Bioinformatics, Volume 39, Issue Supplement_1, June 2023, Pages i185–i193, https://doi.org/10.1093/bioinformatics/btad221

### Example of usage

We obtained the species tree from muti-copy gene family trees using ASTRAL-Pro3 v1.22.3.6 [1] by optimizing the objective function of ASTRAL-Pro [2].
Branch lengths are computed using integrated CASTLES-Pro [3].


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
  - If you see "error" when running `make`, please try `make astral-pro3` instead and file a bug report.
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
* The input trees can have missing taxa, polytomies (unresolved branches), and multi-copy genes.
* When multiple genes from the same species are available, you can ask ASTRAL to force them to be together in the species tree. You can do this in two ways.
  1. You can give multiple genes from the same species the same name in the input gene trees (e.g., `((species_name_A,species_name_B),(species_name_A,species_name_C));`).
  2. OR, a mapping file needs to be provided using the `-a` option. This mapping file should have one line per genes, and each line needs to be in the following formats (e.g., for gene trees like `((gene_A1,gene_B1),(gene_A2,gene_C1));`):

```
gene_A1 species_name_A
gene_A2 species_name_A
gene_B1 species_name_B
gene_B2 species_name_B
gene_B3 species_name_B
...
```

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* (NEW) branch lengths in ***substitution-per-site*** units (IQ-TREE like) for ***all*** branches
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.


# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER and issue the following command:

```
bin/astral-pro3
```

This will give you a list of options available. If you are using Windows, please replace `bin/astral-pro3` with `.\exe\astral-pro3.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
bin/astral-pro3 INPUT_FILE
```
or
```
bin/astral-pro3 -i INPUT_FILE
```

In the first case, INPUT_FILE is ***hard-coded*** to be the ***last argument*** for backward compatibility. 

For example if you want to run `astral-pro3` with input `example/multitree.nw`, then run

```
bin/astral-pro3 example/multitree.nw
```
or
```
bin/astral-pro3 -i example/multitree.nw
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
bin/astral-pro3 -o OUTPUT_FILE INPUT_FILE
```
or
```
bin/astral-pro3 -i INPUT_FILE -o OUTPUT_FILE
```

With `-i INPUT_FILE` option, the order does not matter anymore. For brevity, from here on we will not demonstrate `-i INPUT_FILE` cases.

To save the logs (**also recommended**), run:

```
bin/astral-pro3 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

For example, you can run

```
bin/astral-pro3 -o example/multitree.nw.stree example/multitree.nw 2>example/multitree.nw.log
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
bin/astral-pro3 -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

ASTER also allows rooting at an given outgroup:

```
bin/astral-pro3 --root YOUR_OUTGROUP INPUT_FILE
```

For ASTRAL-Pro, correct rooting is **strongly recommended** to accurately compute branch lengths.

By default, ASTRAL-Pro assumes multiple genes from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/astral-pro3 -a example/multitree_genename.map example/multitree_genename.nw
```

## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
bin/astral-pro3 -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to run with more rounds of placement for ensured optimality, then you can run with
```
bin/astral-pro3 -r 16 -s 16 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```
or simply
```
bin/astral-pro3 -R -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/astral-pro3 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`. If want to score a species tree or you want to place only ***one*** taxon onto the tree, you can use

```
bin/astral-pro3 -r 1 -s 0 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```
or simply,
```
bin/astral-pro3 -C -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/astral-pro3 -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```

Add `-u 0` before `INPUT_FILE` if you want to compute species tree topology only; Add `-u 2` before `INPUT_FILE` if you support and local-PP for all three resolutions of each branch.

```
bin/astral-pro3 -u 0 -o OUTPUT_FILE INPUT_FILE
bin/astral-pro3 -u 2 -o OUTPUT_FILE INPUT_FILE
```

Species tree with more than **5000** taxa may cause **overflow**. Use the following command instead:

```
make astral-pro_int128
bin/astral-pro3_int128 -o OUTPUT_FILE INPUT_FILE
```

If you do not want to compute optimal species tree but instead just want to root and tag gene trees, you can use the following command:

```
bin/astral-pro3 -T -o OUTPUT_FILE INPUT_FILE
```
