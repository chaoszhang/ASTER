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
 
## Bug Reports

Contact ``chaozhang@berkeley.edu``, [``aster-users@googlegroups.com``](https://groups.google.com/forum/#!forum/aster-users), or post on [ASTER issues page](https://github.com/chaoszhang/ASTER/issues).

# Documentations
- The rest of this TUTORIAL file
- [README/astral-pro.md](README/astral-pro.md) for ASTRAL and ASTRAL-Pro; [README/wastral.md](README/wastral.md) for weighted ASTRAL series; [README/asterisk.md](README/asterisk.md) for ASTERISK series
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
1. [ASTRAL](tutorial/astral.md)
2. [ASTRAL-Pro](tutorial/astral-pro.md)
3. [Weighted ASTRAL by Branch Support](README/wastral.md)
4. [Weighted ASTRAL by Branch Length](README/wastral.md)
5. [Weighted ASTRAL - Hybrid](README/wastral.md)
6. [ASTERISK](README/asterisk.md)

# ACKNOWLEGEMENT
ASTER code uses Regularized Incomplete Beta Function by Lewis Van Winkle under zlib License. Code is contributed by Chao Zhang supervised by Siavash Mirarab.
