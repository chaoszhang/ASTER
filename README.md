# Accurate Species Tree EstimatoR (ASTER❋)

[<img src="misc/ASTER.png" width="250"/>](misc/ASTER.png)

A family of optimatization algorithms for species tree inference:
1. [ASTRAL-IV](tutorial/astral4.md) (from unrooted gene tree ***topologies*** with integrated CASTLES-II for branch lengths)
2. [ASTRAL-Pro3](tutorial/astral-pro3.md) (from unrooted gene ***family*** tree topologies with integrated CASTLES-Pro)
3. [Weighted ASTRAL](tutorial/wastral.md) (from unrooted gene trees with branch ***lengths*** and/or ***supports***)
4. [CASTER-site](tutorial/caster-site.md) (from ***whole genome alignments*** or aligned sequences)
5. [CASTER-pair](tutorial/caster-site.md) (from ***whole genome alignments*** or aligned sequences)
6. [WASTER](tutorial/waster.md) (from ***raw reads***)
7. SISTER (from optical-map-like distance data or shape data)
8. MONSTER

# Announcements

## Available on Bioconda

We thank Dirk-Jan M. van Workum, Joshua Zhuang, and many other contributors who make ASTER available on [Bioconda](https://bioconda.github.io/recipes/aster/README.html).

## Integrated in Phylosuite

Many ASTER tools have been integrated in [PhyloSuite](http://phylosuite.jushengwu.com/), an integrated and scalable desktop platform for streamlined molecular sequence data management and evolutionary phylogenetics studies.

## GUI for Windows users

Please check out our software with GUI for ASTRAL, ASTRAL-Pro, and wASTRAL. Simply download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip), extract the contents, enter `exe` folder, and click `aster-gui.exe`.
 
## Bug Reports

Contact ``chao.zhang@sund.ku.dk``, [``aster-users@googlegroups.com``](https://groups.google.com/forum/#!forum/aster-users), or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this README file
- Program specific tutorials (see EXECUTION section)
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
    - Unix (MacOS) users should be prompted for installing `g++` and please click "install". If no prompt, try `g++`
2. Binary files should be in the `bin` folder.

## For Windows users
- Executables for x86-64 are available in `exe` folder and it is **very likely** that they already work.
- [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) is HIGHLY recommanded if you need to install on your own! Please follow instructions in "For Linux/Unix/WSL users" section.
- To compile windows excutables:
  1. Download [MinGW](https://sourceforge.net/projects/mingw-w64/) and install ***posix*** version for your architecture (eg. x86-64)
  2. Add path to `bin` folder of MinGW to [system environment variable `PATH`](https://www.google.com/search?q=Edit+the+system+environment+variables+windows)
  3. Double click `make.bat` inside the downloaded directory

# EXECUTION
Please click the link below:
1. [ASTRAL-IV](tutorial/astral4.md)
2. [ASTRAL-Pro3](tutorial/astral-pro3.md)
3. [Weighted ASTRAL](tutorial/wastral.md)
4. [CASTER-site](tutorial/caster-site.md)
5. [CASTER-pair](tutorial/caster-site.md)
6. [WASTER](tutorial/waster.md)

# HELP ME CHOOSE A SUITABLE TOOL

Q: I want to reconstruct a phylogenetic tree without alignments.

A: I recommend **WASTER**, which can accurately infer family-level species tree with short reads with even 1.5X coverage. This tool can also be used to build an adequate order-level guide tree for CACTUS. (You will be amazed by how accurate this "guide tree" is!)

Q: I have a supermatrix of SNPs in fasta/phylip format and I want a "quick-and-dirty" run to get an adequate phylogenetic tree.

A: I recommend **CASTER-site**, which is usually 1-2 magnitudes faster than concatenation-based maximum likelihood methods yet more accurate in presence of incomplete lineage sorting with enough data.

Q: My dataset has a lot of muti-copy genes (e.g. plants) and I want to make an effort to utilize these precious signals.

A: I highly recommend **ASTRAL-Pro3**, which takes as input non-rooted non-labelled gene family trees. ASTRAL-Pro3 does not need to know the homology relationships of genes, but you still need to reconstruct gene family trees by yourself using RAxML/IQTree/Fasttree.

Q: I have aligned genomes (>10M sites) and the average nucleotide identity is >80% between closely related species (e.g. birds, mammals, or abundant taxon sampling).

A: I recommend **CASTER-site** (faster) and **CASTER-pair** (slower). Those methods are usually 1-2 magnitudes faster than concatenation-based maximum likelihood methods yet more accurate in presence of incomplete lineage sorting. Please run both and select the species tree that makes more sense.

Q: I have gene trees with branch lengths and Bootstrap/Baysian supports and I know that horizontal gene transfers and hybridizations are rare.

A: I recommend **Weighted ASTRAL**. It utilizes branch lengths and supports to improve accuracy.

Q: I have gene trees but they do not satisfy the requirements for wASTRAL.

A: You can still use **ASTRAL-IV**. By the way, ASTRAL-IV is also useful for finding the supertree.

# ACKNOWLEGEMENT
ASTER code uses Regularized Incomplete Beta Function by Lewis Van Winkle under zlib License. Code is contributed by Chao Zhang supervised by Siavash Mirarab and Rasmus Nielsen.
