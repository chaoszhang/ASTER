# Accurate Species Tree EstimatoR
A family of Optimatization algorithms for species tree inference:
1. [ASTRAL](README/astral-pro.md) (re-implemented in C++)
2. [ASTRAL-Pro](README/astral-pro.md) (re-implemented in C++ for better running time, memory consumption, and accuracy)
3. [Weighted ASTRAL by Branch Support](README/wastral.md)
4. [Weighted ASTRAL by Branch Length](README/wastral.md)
5. [Weighted ASTRAL - Hybrid](README/wastral.md)
6. [ASTERISK](README/asterisk.md)

# Bug Reports:

Contact ``aster-users@googlegroups.com`` or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this README file
- [README/astral-pro.md](README/astral-pro.md) for ASTRAL and ASTRAL-Pro; [README/wastral.md](README/wastral.md) for weighted ASTRAL series; [README/asterisk.md](README/asterisk.md) for ASTERISK series
- Forums:
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER issues page](https://github.com/chaoszhang/ASTER/issues)

# INSTALLATION
For most users, installing ASTER is VERY easy!
Download using one of two approaches:
  - You simply need to download the [zip file](https://github.com/chaoszhang/ASTER/archive/refs/heads/master.zip) and extract the contents to a folder of your choice.
  - Alternatively, you can clone the [github repository](https://github.com/chaoszhang/ASTER.git).

## For Linux/Unix/WSL users
1. In terminal, `cd` into the downloaded directory and run `make`.
  - If you see `*** Installation complete! ***` then you are done!
  - If you see `Command 'g++' not found` then
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
- Executables for x86-64 are available in `exe` folder
- [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install) is HIGHLY recommanded if you need to install on your own! Please follow instructions in "For Linux/Unix/WSL users" section.
- To compile windows excutables:
  1. Download [MinGW](https://sourceforge.net/projects/mingw-w64/) and install ***posix*** version for your architecture (eg. x86-64)
  2. Add path to `bin` folder of MinGW to [system environment variable `PATH`](https://www.google.com/search?q=Edit+the+system+environment+variables+windows)
  3. Double click `make.bat` inside the downloaded directory

# EXECUTION
ASTER currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to [`bin`](Linux/Unix/WSL) or [`exe`](Win) in the location where you have downloaded the software, find the name of your program of interest, and issue the following command:

```
./PROGRAM_NAME 
```

This will give you a list of options available. Windows program name should include `.exe` extension.

To find the species tree with input from in a file called `INPUT_FILE`, use:

```
./PROGRAM_NAME INPUT_FILE
```

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
