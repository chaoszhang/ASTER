# Accurate Species Tree EstimatoR (ASTER❋)

[<img src="logo.png" width="250"/>](logo.png)

A family of optimatization algorithms for species tree inference:
1. [ASTRAL](tutorial/astral.md) (re-implemented in C++, suitable for large data, multi-individual, and super-tree)
2. [ASTRAL-Pro](tutorial/astral-pro.md) (re-implemented in C++ for better running time, memory consumption, and usability)
3. [Weighted ASTRAL by Branch Support](README/wastral.md)
4. [Weighted ASTRAL by Branch Length](README/wastral.md)
5. [Weighted ASTRAL - Hybrid](README/wastral.md)
6. [ASTERISK](README/asterisk.md)

## GUI for Windows users (NEW)

Please check out our software with GUI. Simply download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/Windows.zip), extract the contents, enter `exe` folder, and click `aster-gui.exe`. 
<!--
GUI for other plantforms is also available at `src/aster-gui.py` with dependencies:
```
pip3 install pyqt5 pyqt5-tools
```
-->
 
## Bug Reports:

Contact ``aster-users@googlegroups.com`` or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this README file
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

# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to `bin`(Linux/Unix/WSL) or `exe`(Win) in the location where you have downloaded the software, find the name of your program of interest, and issue the following command:

```
./PROGRAM_NAME 
```

This will give you a list of options available. Windows program name should include `.exe` extension.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
./PROGRAM_NAME INPUT_FILE
```

Currently, INPUT_FILE is hard-coded to be the ***last argument***. 

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option before `INPUT_FILE`(**Strongly recommended**):

```
./PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE
```

To save the logs (**also recommended**), run:

```
./PROGRAM_NAME -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

ASTER supports multi-threading. To run program with 4 threads, add `-t 4` before `INPUT_FILE`:

```
./PROGRAM_NAME -t 4 -o OUTPUT_FILE INPUT_FILE 2>LOG_FILE
```

## Advanced Options

ASTER algorithm first performs `R` (4 by default) rounds of search and then repeatedly performs `S` (4 by default) rounds of subsampling and exploration until no improvement found.

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

Add `-u 0` before `INPUT_FILE` if you want to compute species tree topology only; Add `-u 2` before `INPUT_FILE` if you support and local-PP for all three resolutions of each branch for **supported programs**.

```
./PROGRAM_NAME -u 0 -o OUTPUT_FILE INPUT_FILE
./PROGRAM_NAME -u 2 -o OUTPUT_FILE INPUT_FILE
```
# INPUT, OUTPUT, AND PROGRAM SPECIFIC FEATURES
See 
 - [README/astral-pro.md](README/astral-pro.md) for ASTRAL and ASTRAL-Pro;
 - [README/wastral.md](README/wastral.md) for weighted ASTRAL series;
 - [README/asterisk.md](README/asterisk.md) for ASTERISK series.

# ACKNOWLEGEMENT
ASTER code uses Regularized Incomplete Beta Function by Lewis Van Winkle under zlib License. Code is contributed by Chao Zhang supervised by Siavash Mirarab.
