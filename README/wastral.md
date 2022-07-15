# Weighted ASTRAL Series (wASTRAL)
1. Weighted ASTRAL by Branch Support (astral-weighted)
2. Weighted ASTRAL by Branch Length (astral-lengthweighted)
3. Weighted ASTRAL - Hybrid (astral-hybrid)

## Publication

TBD

## Bug Reports

Contact ``aster-users@googlegroups.com`` or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this file
- [README.md](../README.md)
- Forums:
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER issues page](https://github.com/chaoszhang/ASTER/issues)

# INSTALLATION
See [README.md](../README.md).

# EXECUTION
Weighted ASTRAL currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to `bin`(Linux/Unix/WSL) or `exe`(Win) in the location where you have downloaded the software, find the name of your program of interest, and issue the following command:

```
./PROGRAM_NAME
```

Replace `PROGRAM_NAME` with `astral-weighted` (support), `astral-lengthweighted` (length), or `astral-hybrid` (hybrid). This will give you a list of options available. On Windows, replace `./PROGRAM_NAME` with `.\PROGRAM_NAME.exe`.

To find the species tree with input gene trees from in a file called `INPUT_FILE`, use:

```
./PROGRAM_NAME INPUT_FILE
```

Currently, INPUT_FILE is hard-coded to be the ***last argument***. 

***Notice:*** For `astral-weighted` (support) and `astral-hybrid` (hybrid), you ***must*** provide max (`-x`) and min (`-n`) of support value to the program; otherwise it will default to `max=100, min=0`. For **Bootstrap** support, the default value should be fine; for **local Baysian** support, `-x 1 -n 0.333` is recommended; for other **probability & likelihood** support, `-x 1 -n 0` may be more reasonable. Please ensure that you provide the correct max and min; otherwise, the output of the program ***may not be meaningful***.

```
./PROGRAM_NAME -x MAX_SUPPORT -n MIN_SUPPORT INPUT_FILE
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
./PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE
```

To save the logs (**also recommended**), run:

```
./PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

Weighted ASTRAL supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
./PROGRAM_NAME -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

Example: 
```
./astral-hybrid -t 4 -o ../example/genetree.astral-hybrid.nw ../example/genetree.nw 2>../example/genetree.astral-hybrid.log
```

## Advanced Options

Weighted ASTRAL algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
./PROGRAM_NAME -t T -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
``` 

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
./PROGRAM_NAME -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`.

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
./PROGRAM_NAME -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```

If you want **long double** instead **double** precision. Use the following command instead (***It is slow!***):

```
./PROGRAM_NAME_precise -o OUTPUT_FILE INPUT_FILE
```

**Currently branch lengths and local-PPs are only for `astral-hybrid`.** Add `-u 0` before `INPUT_FILE` if you want to compute species tree topology only; Add `-u 2` before `INPUT_FILE` if you support and local-PP for all three resolutions of each branch.

```
./PROGRAM_NAME -u 0 -o OUTPUT_FILE INPUT_FILE
./PROGRAM_NAME -u 2 -o OUTPUT_FILE INPUT_FILE
```

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

1. **Weighted ASTRAL by Branch Support (astral-weighted)**: Non-root interal node labels must be a non-negative number being support. Eg. `((A,B)100,(C,D)0);` or `((A:1,B:1)1.0:1,(C:1,D:1)0.333:0);`.
2. **Weighted ASTRAL by Branch Length (astral-lengthweighted)**: Non-root labels must have branch lengths after `:`. Eg. `((A,B):1,(C,D):0);` or `((A:1,B:1)1.0:1,(C:1,D:1)0.333:0);`.
3. **Weighted ASTRAL - Hybrid (astral-hybrid)**: Non-root interal node labels must be a non-negative number being support before `:` and non-root labels must have branch lengths after `:`. Eg. `((A:1,B:1)100:1,(C:1,D:1)0:0);` or `((A:1,B:1)1.0:1,(C:1,D:1)0.333:0);`.

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch lengths in coalescent units for astral-weighted and in **combined (coalensecent + 2 * substitution) units** (only for internal branches)
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.

The weighted ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).
