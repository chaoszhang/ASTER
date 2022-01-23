mkdir exe
g++ -std=gnu++11 -O2 -pthread src/astral.cpp -o exe/astral.exe
g++ -std=gnu++11 -O2 -pthread src/astral-pro.cpp -o exe/astral-pro.exe
g++ -std=gnu++11 -O2 -pthread src/astral-weighted.cpp -o exe/astral-weighted.exe
g++ -std=gnu++11 -O2 -pthread src/astral-lengthweighted.cpp -o exe/astral-lengthweighted
g++ -std=gnu++11 -O2 -pthread src/astral-hybrid.cpp -o exe/astral-hybrid.exe
g++ -std=gnu++11 -O2 -pthread src/asterisk.cpp -o exe/asterisk.exe