#ifndef TUTORIAL
#define TUTORIAL

#include<regex>
#include<fstream>

class MDGenerator{
const string APRO_UNIQUE_INTRO = R"V0G0N(# Accurate Species Tree ALgorithm for PaRalogs and Orthologs (ASTRAL-Pro3)
ASTRAL-Pro stands for ASTRAL for PaRalogs and Orthologs. ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees and is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL-pro extends ASTRAL to allow multi-copy genes. ASTRAL-pro finds the species tree that has the maximum number of shared induced quartet tree equivalent classes with the set of gene trees, subject to the constraint that the set of bipartitions in the species tree comes from a predefined set of bipartitions. Please see the paper below for the definition of the PL-quartet scores, which is what ASTRAL-Pro optimizes. We refer to the tool both as A-Pro and ASTRAL-Pro. 

ASTRAL-Pro3 re-implements [ASTRAL-Pro](https://github.com/chaoszhang/A-pro) in an equally accurate yet **faster**, and **easier to install** and **lower memory consumption** way.
ASTRAL-Pro3 also integrates [CASTLES-Pro](https://github.com/ytabatabaee/CASTLES) and thus computes terminal and internal branch lengths in substitution-per-site units.

## Publication

[1] Chao Zhang, Siavash Mirarab, ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees, Bioinformatics, 2022, btac620, https://doi.org/10.1093/bioinformatics/btac620

[2] Chao Zhang, Celine Scornavacca, Erin K Molloy, Siavash Mirarab, ASTRAL-Pro: Quartet-Based Species-Tree Inference despite Paralogy, Molecular Biology and Evolution, Volume 37, Issue 11, November 2020, Pages 3292–3307, https://doi.org/10.1093/molbev/msaa139

[3] Yasamin Tabatabaee, Chao Zhang, Tandy Warnow, Siavash Mirarab, Phylogenomic branch length estimation using quartets, Bioinformatics, Volume 39, Issue Supplement_1, June 2023, Pages i185–i193, https://doi.org/10.1093/bioinformatics/btad221

### Example of usage

We obtained the species tree from muti-copy gene family trees using ASTRAL-Pro3 VERSION [1] by optimizing the objective function of ASTRAL-Pro [2].
Branch lengths are computed using integrated CASTLES-Pro [3].

)V0G0N";

const string ASTRAL_UNIQUE_INTRO = R"V0G0N(# Accurate Species Tree ALgorithm (ASTRAL-IV)
ASTRAL is a tool for estimating an unrooted species tree given a set of unrooted gene trees. ASTRAL is statistically consistent under the multi-species coalescent model (and thus is useful for handling incomplete lineage sorting, i.e., ILS). ASTRAL finds the species tree that has the maximum number of shared induced quartet trees with the set of gene trees, subject to the constraint that the set of tripartitions in the species tree comes from a predefined set of tripartitions.

ASTRAL-IV re-implements [ASTRAL](https://github.com/smirarab/ASTRAL) as a scalable alternative to ASTRAL on datasets for which ASTRAL is not suitable (e.g. large datasets, multi-individual, and gene trees with missing taxa).
ASTRAL-IV also integrates [CASTLES-II](https://github.com/ytabatabaee/CASTLES) and thus computes terminal and internal branch lengths in substitution-per-site units.

As a scalable alternative to ASTRAL-III, ASTRAL-IV lacks of some features of ASTRAL-III (e.g. bootstrapping). You can work around by first computing optimal tree with ASTRAL-IV and use the ASTRAL-IV output tree as `-q` option to ASTRAL-III. 

## Publication

[1] Chao Zhang, Siavash Mirarab, Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees, Molecular Biology and Evolution, 2022, msac215, https://doi.org/10.1093/molbev/msac215

[2] Chao Zhang, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).

[3] Yasamin Tabatabaee, Chao Zhang, Tandy Warnow, Siavash Mirarab, Phylogenomic branch length estimation using quartets, Bioinformatics, Volume 39, Issue Supplement_1, June 2023, Pages i185–i193, https://doi.org/10.1093/bioinformatics/btad221

### Example of usage

We obtained the species tree from gene trees using ASTRAL-IV VERSION [1] by optimizing the objective function of ASTRAL [2].
Branch lengths are computed using integrated CASTLES-II [3].

)V0G0N";

const string WASTRAL_UNIQUE_INTRO = R"V0G0N(# Weighted ASTRAL (wASTRAL)
Weighted ASTRAL program has the following modes for different weighting schemes:
1. Weighted ASTRAL by Branch Support (mode 2)
2. Weighted ASTRAL by Branch Length (mode 3)
3. Weighted ASTRAL - Hybrid (default)

Weighted ASTRAL series introduce threshold-free weighting schemes for the quartet-based species tree inference, the metric used in the popular method ASTRAL. By reducing the impact of quartets with low support or long terminal branches (or both), weighting provides stronger theoretical guarantees and better empirical performance than the unweighted ASTRAL. Our results show that weighted ASTRAL improves the utility of summary methods and can reduce the incongruence often observed across analytical pipelines.

## Publication

[1] Chao Zhang, Siavash Mirarab, Weighting by Gene Tree Uncertainty Improves Accuracy of Quartet-based Species Trees, Molecular Biology and Evolution, 2022, msac215, https://doi.org/10.1093/molbev/msac215

(OPTIONAL) [2] Chao Zhang, Maryam Rabiee, Erfan Sayyari, and Siavash Mirarab. 2018. “ASTRAL-III: Polynomial Time Species Tree Reconstruction from Partially Resolved Gene Trees.” BMC Bioinformatics 19 (S6): 153. [doi:10.1186/s12859-018-2129-y](https://doi.org/10.1186/s12859-018-2129-y).

### Example of usage

We obtained the species tree from gene trees using wASTRAL VERSION [1]. 
(OPTIONAL) This method reimplements ASTRAL-III [2] and takes into account phylogenetic uncertainty by intergrating signals from branch length and branch support in gene trees.

)V0G0N";

const string CASTER_UNIQUE_INTRO = R"V0G0N(# Coalescence-aware Alignment-based Species Tree EstimatoR (CASTER)

[<img src="../misc/CASTER.png" width="500"/>](../misc/CASTER.png)

Genome-wide data have the promise of dramatically improving phylogenetic inferences. Yet, inferring the true phylogeny remains a challenge, mainly because the evolutionary histories of different genomic regions differ. The traditional concatenation approach ignores such differences, resulting in both theoretical and empirical shortcomings. In response, many discordance-aware inference methods have been developed. Yet, all have their own weaknesses. Many methods rely on short recombination-free genomic segments to build gene trees and thus suffer from a lack of signals for gene tree reconstruction, resulting in poor species tree. Some methods wrongly assume that the rate of evolution is uniform across the species tree. Yet, others lack enough scalability to analyze phylogenomic data.

We introduce a new site-based species tree inference method that seeks to address these challenges without reconstructing gene trees. Our method, called CASTER (Coalescence-aware Alignment-based Species Tree EstimatoR), has two flavors: CASTER-site and CASTER-pair. The first is based on patterns in individual sites and the second is based on pairs of sites.

CASTER has several outstanding features:
1. CASTER introduces two new optimization objectives based on genomic site patterns of four species; we show that optimizing these objectives produces two estimators: CASTER-site is statistically consistent under MSC+HKY model while allowing mutation rate to change across sites and across species tree branches; CASTER-pair is statistically consistent under MSC+GTR model under further assumptions.
2. CASTER comes with a scalable algorithm to optimize the objectives summed over all species quartets. Remarkably, its time complexity is linear to the number of sites and at most quasi-quadratic with respect to the number of species.
3. CASTER can handle multiple samples per species, and CASTER-site specifically can work with allele frequencies of unphased multiploid.
4. CASTER is extremenly memory efficent, requiring <1 byte per SNP per sample

Under extensive simulation of genome-wide data, including recombination, we show that both CASTER-site and CASTER-pair out-perform concatenation using RAxML-ng, as well as discordance-aware methods SVDQuartets and wASTRAL in terms of both accuracy and running time. Noticeably, CASTER-site is 60–150X faster than the alternative methods. It reconstructs an Avian tree of 51 species from aligned genomes with 254 million SNPs in only 3.5 hours on an 8-core desktop machine with 32 GB memory. It can also reconstruct a species tree of 201 species with approximately 2 billion SNPs using a server of 256 GB memory.

Our results suggest that CASTER-site and CASTER-pair can fulfill the need for large-scale phylogenomic inferences.

## Publication

Chao Zhang, Rasmus Nielsen, Siavash Mirarab, CASTER: Direct species tree inference from whole-genome alignments. Science (2025) https://www.science.org/doi/10.1126/science.adk9688

## Notice

Since CASTER-site and CASTER-pair assume different models, please run both and choose the result that makes more sense if you can.

)V0G0N";

const string WASTER_UNIQUE_INTRO = R"V0G0N(# Without-Alignment/Assembly Species Tree EstimatoR † (WASTER)

[<img src="../misc/WASTER.png" width="500"/>](../misc/WASTER.png)

WASTER is a coalesence-aware ***de novo*** species tree inference tool, which means it can take as inputs raw reads in FASTQ format.
Paticularly, WASTER can accurately infer species tree even from Illumina reads with only ***1.5X depth***.
WASTER infers the species tree by first calling SNPs from reads/assembies and then invoking CASTER to reconstruct the species tree from the SNPs.

## Publication

Chao Zhang, Rasmus Nielsen, WASTER: Practical de novo phylogenomics from low-coverage short reads, bioRxiv (2025) https://doi.org/10.1101/2025.01.20.633983

)V0G0N";

const string SHARED_INTRO = R"V0G0N(
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
  - If you see "error" when running `make`, please try `make PROGRAM_NAME` instead and file a bug report.
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
)V0G0N";

const string APRO_IO = R"V0G0N(
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

)V0G0N";
// The ASTRAL-Pro tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).

const string ASTRAL_IO = R"V0G0N(
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
* (NEW) branch lengths in ***substitution-per-site*** units (IQ-TREE like) for ***all*** branches
* branch supports measured as [local posterior probabilities](http://mbe.oxfordjournals.org/content/early/2016/05/12/molbev.msw079.short?rss=1)
* It can also annotate branches with other quantities, such as quartet supports and localPPs for all three topologies.

)V0G0N";
// The ASTRAL tree leaves the branch length of terminal branches empty. Some tools for visualization and tree editing do not like this (e.g., ape). In FigTree, if you open the tree several times, it eventually opens up (at least on our machines). In ape, if you ask it to ignore branch lengths all together, it works. In general, if your tool does not like the lack of terminal branches, you can add a dummy branch length, [as in this script](https://github.com/smirarab/global/blob/master/src/mirphyl/utils/add-bl.py).

const string WASTRAL_IO = R"V0G0N(
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
)V0G0N";

const string CASTER_IO = R"V0G0N(
# STOP!
Please make sure you removed paralogous alignment regions using `RepeatMasker` or alike. This will improve the accuracy of CASTER on biological datasets!

A detailed [walkthrough](../misc/caster-linux-walkthrough.md) is also available (English/中文).

A walkthrough for [sliding window](../misc/caster-sliding-window-walkthrough.md) analysis is also available.

# INPUT
* The input by default is a single MSA in Fasta format
* The input can also be a text file containing a list of Fasta files (one file per line) if you add `-f list` to input arguments
* The input can also be a single Phylip files or vertically concatenated Phylip files in one file by adding `-f phylip` to input arguments
* The input files can have missing taxa and multiple individuals/copies per species.
* When individuals/copies/genes from the same species are available, you need to let CASTER to know that they are from the same species. You can do this in two ways.
  1. You can give multiple individuals/copies/genes from the same species the same name in the input alignment.
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

Examples:

Single Fasta file:
```
>species_A
AAA
>species_C
CCC
>species_G
GGG
>species_T
TTT
```

Single Phylip file, multiple individuals with mapping file:
```
5 3
individual_A1 AAA
individual_A2 AAA
species_C CCC
species_G GGG
species_T TTT
```
Mapping file:
```
individual_A1 species_A
individual_A2 species_A
```

Multiple Fasta file:
```
gene1.fasta
gene2.fasta
```
In `gene1.fasta`:
```
>species_A
AAA
>species_C
CCC
>species_G
GGG
>species_T
TTT
```
In `gene2.fasta` (order can change, fasta files can have missing taxa):
```
>species_C
CC
>species_A
AA
>species_T
TT
```

Multiple Phylip file, multiploid without mapping file (if using `CASTER-pair`, genes must be phased; if using `CASTER-site`, you can arbitrarily phase them):
```
6 3
species_A AAA
species_A AAA
species_A AAA
species_A AAA
species_C CCC
species_C CCC
4 2
species_A AA
species_A AA
species_T TT
species_T TT
```

Notice: only `CASTER-site` works on unphased SNPs, you can translate VCF files into Fasta (or Phylip) in the following way.

VCF:
```
species_A
A/A/C/T
A/A/G/G
```
Fasta (order and phasing do not matter):
```
>species_A
AA
>species_A
AA
>species_A
CG
>species_A
TG
```

# OUTPUT
The output in is Newick format and gives:

* the species tree topology
* branch supports measured as local bootstrap support (>95.0 means good)
* It can also annotate branches with other quantities, such as quartet scores and local bootstraps for all three topologies.
)V0G0N";

const string WASTER_IO = R"V0G0N(
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
)V0G0N";

const string SHARED_EXE = R"V0G0N(
# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER and issue the following command:

```
bin/PROGRAM_NAME
```

This will give you a list of options available. If you are using Windows, please replace `bin/PROGRAM_NAME` with `.\exe\PROGRAM_NAME.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
bin/PROGRAM_NAME INPUT_FILE
```
or
```
bin/PROGRAM_NAME -i INPUT_FILE
```

In the first case, INPUT_FILE is ***hard-coded*** to be the ***last argument*** for backward compatibility. 

For example if you want to run `PROGRAM_NAME` with input `example/EXAMPLE_INPUT`, then run

```
bin/PROGRAM_NAME example/EXAMPLE_INPUT
```
or
```
bin/PROGRAM_NAME -i example/EXAMPLE_INPUT
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
bin/PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE
```
or
```
bin/PROGRAM_NAME -i INPUT_FILE -o OUTPUT_FILE
```

With `-i INPUT_FILE` option, the order does not matter anymore. For brevity, from here on we will not demonstrate `-i INPUT_FILE` cases.

To save the logs (**also recommended**), run:

```
bin/PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

For example, you can run

```
bin/PROGRAM_NAME -o example/EXAMPLE_INPUT.stree example/EXAMPLE_INPUT 2>example/EXAMPLE_INPUT.log
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
bin/PROGRAM_NAME -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

ASTER also allows rooting at an given outgroup:

```
bin/PROGRAM_NAME --root YOUR_OUTGROUP INPUT_FILE
```
)V0G0N";

const string APRO_UNIQUE_EXE = R"V0G0N(
For ASTRAL-Pro, correct rooting is **strongly recommended** to accurately compute branch lengths.

By default, ASTRAL-Pro assumes multiple genes from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/astral-pro3 -a example/multitree_genename.map example/multitree_genename.nw
```
)V0G0N";

const string ASTRAL_UNIQUE_EXE = R"V0G0N(
For ASTRAL, correct rooting is **strongly recommended** to accurately compute branch lengths.

By default, ASTRAL assumes multiple individuals from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/astral4 -a example/genetree.map example/genetree.nw
```

When your dataset has no more than 50 species and no more than 500 genes, you may want to run with more rounds using `-R` (see below). 
)V0G0N";

const string WASTRAL_UNIQUE_EXE = R"V0G0N(
By default, wASTRAL assumes multiple individuals/alleles from the same species in the same input gene trees having the same name. Alternatively, a mapping file needs to be provided using the `-a` option (see INPUT section). For example,

```
bin/wastral -a example/genetree.map example/genetree.nw
```

When your dataset has no more than 50 species and no more than 500 genes, you may want to run with more rounds using `-R` (see below). 
)V0G0N";

const string WASTER_UNIQUE_EXE = R"V0G0N(
***Notice:*** By default, WASTER will create a 32 GB look-up table to find SNPs. So, if you are running on a machine with <64 GB memory, you need to shrink the look-up table size using `-k 8` or `-k 7` even for the example run.
Using `-k 8` requires about 4 GB memory and `-k 7` only requires <1 GB memory.
Try the following example run:
```
bin/waster-site -k 7 example/waster/input_list.txt
```
)V0G0N";

const string SHARED_ADV = R"V0G0N(
## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

```
bin/PROGRAM_NAME -r R -s S -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to run with more rounds of placement for ensured optimality, then you can run with
```
bin/PROGRAM_NAME -r 16 -s 16 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```
or simply
```
bin/PROGRAM_NAME -R -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

If you want to place taxa on an existing ***fully resolved*** species tree, you can use `-c SPECIES_TREE_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/PROGRAM_NAME -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

Specifically, you can score and annotate a ***fully resolved*** species tree containing all taxa with `-c SPECIES_TREE_IN_NEWICK_FORMAT`. If want to score a species tree or you want to place only ***one*** taxon onto the tree, you can use

```
bin/PROGRAM_NAME -r 1 -s 0 -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```
or simply,
```
bin/PROGRAM_NAME -C -o OUTPUT_FILE -c SPECIES_TREE_IN_NEWICK_FORMAT INPUT_FILE
```

If you want to give hints by providing candidate species trees or trees similar to the species tree, you can use `-g SPECIES_TREES_IN_NEWICK_FORMAT` before `INPUT_FILE`:

```
bin/PROGRAM_NAME -o OUTPUT_FILE -g SPECIES_TREES_IN_NEWICK_FORMAT INPUT_FILE
```
)V0G0N";

const string APRO_UNIQUE_ADV = R"V0G0N(
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
)V0G0N";

const string ASTRAL_UNIQUE_ADV = R"V0G0N(
Add `-u 0` before `INPUT_FILE` if you want to compute species tree topology only; Add `-u 2` before `INPUT_FILE` if you support and local-PP for all three resolutions of each branch.

```
bin/astral4 -u 0 -o OUTPUT_FILE INPUT_FILE
bin/astral4 -u 2 -o OUTPUT_FILE INPUT_FILE
```

Species tree with more than **5000** taxa may cause **overflow**. Use the following command instead:

```
make astral_int128
bin/astral4_int128 -o OUTPUT_FILE INPUT_FILE
```
)V0G0N";

const string WASTRAL_UNIQUE_ADV = R"V0G0N(
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

)V0G0N";
    const string APRO = "astral-pro3";
    const string ASTRAL = "astral4";
    const string WASTRAL = "wastral";
    const string CASTER = "caster-site";
    const string WASTER = "waster";

    string replace(string txt, string from, string to){
        return regex_replace(txt, regex(from), to);
    }

    string uniqueIntro(){
        string rawtext = "# Accurate Species Tree EstimatoR (ASTER❋)";
        if (programName == APRO) rawtext = APRO_UNIQUE_INTRO;
        if (programName == ASTRAL) rawtext = ASTRAL_UNIQUE_INTRO;
        if (programName == WASTRAL) rawtext = WASTRAL_UNIQUE_INTRO;
        if (programName == CASTER) rawtext = CASTER_UNIQUE_INTRO;
        if (programName == WASTER) rawtext = WASTER_UNIQUE_INTRO;
        return replace(rawtext, "VERSION", version);
    }

    string sharedIntro(){
        if (programName == APRO) return replace(SHARED_INTRO, "PROGRAM_NAME", "astral-pro3");
        if (programName == ASTRAL) return replace(SHARED_INTRO, "PROGRAM_NAME", "astral4");
        if (programName == WASTRAL) return replace(SHARED_INTRO, "PROGRAM_NAME", "wastral");
        if (programName == CASTER) return replace(SHARED_INTRO, "PROGRAM_NAME", "caster-site");
        if (programName == WASTER) return replace(SHARED_INTRO, "PROGRAM_NAME", "waster");
        return SHARED_INTRO;
    }

    string inputOutput(){
        if (programName == APRO) return APRO_IO;
        if (programName == ASTRAL) return ASTRAL_IO;
        if (programName == WASTRAL) return WASTRAL_IO;
        if (programName == CASTER) return CASTER_IO;
        if (programName == WASTER) return WASTER_IO;
        return "";
    }

    string sharedExe(){
        if (programName == APRO) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "astral-pro3"), "EXAMPLE_INPUT", "multitree.nw");
        if (programName == ASTRAL) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "astral4"), "EXAMPLE_INPUT", "genetree.nw");
        if (programName == WASTRAL) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "wastral"), "EXAMPLE_INPUT", "genetree.nw");
        if (programName == CASTER) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "caster-site"), "EXAMPLE_INPUT", "genetrees.tre_1.fas");
        if (programName == WASTER) return replace(replace(SHARED_EXE, "PROGRAM_NAME", "waster"), "EXAMPLE_INPUT", "waster/input_list.txt");
        return SHARED_EXE;
    }

    string uniqueExe(){
        if (programName == APRO) return APRO_UNIQUE_EXE;
        if (programName == ASTRAL) return ASTRAL_UNIQUE_EXE;
        if (programName == WASTRAL) return WASTRAL_UNIQUE_EXE;
        if (programName == WASTER) return WASTER_UNIQUE_EXE;
        return "";
    }

    string sharedAdv(){
        if (programName == APRO) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "astral-pro3"), "EXAMPLE_INPUT", "multitree.nw");
        if (programName == ASTRAL) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "astral4"), "EXAMPLE_INPUT", "genetree.nw");
        if (programName == WASTRAL) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "wastral"), "EXAMPLE_INPUT", "genetree.nw");
        if (programName == CASTER) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "caster-site"), "EXAMPLE_INPUT", "genetrees.tre_1.fas");
        if (programName == WASTER) return replace(replace(SHARED_ADV, "PROGRAM_NAME", "waster"), "EXAMPLE_INPUT", "waster/input_list.txt");
        return SHARED_ADV;
    }

    string uniqueAdv(){
        if (programName == APRO) return APRO_UNIQUE_ADV;
        if (programName == ASTRAL) return ASTRAL_UNIQUE_ADV;
        if (programName == WASTRAL) return WASTRAL_UNIQUE_ADV;
        return "";
    }

    string programName;

public:
    
    static string version;

    MDGenerator(string programName): programName(programName){
        string ui = uniqueIntro();
        string si = sharedIntro();
        string io = inputOutput();
        string se = sharedExe();
        string ue = uniqueExe();
        string sa = sharedAdv();
        string ua = uniqueAdv();
        ofstream fout(string("tutorial/") + programName + ".md");
        fout << (ui + si + io + se + ue + sa + ua);
    }

};

string MDGenerator::version;

#endif
