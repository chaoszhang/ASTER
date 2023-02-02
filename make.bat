mkdir exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral.cpp -o exe/astral.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral.cpp -o exe/astral_int128.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-pro.cpp -o exe/astral-pro.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-pro.cpp -o exe/astral-pro_int128.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-weighted.cpp -o exe/astral-weighted.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-weighted.cpp -o exe/astral-weighted_precise.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-lengthweighted.cpp -o exe/astral-lengthweighted.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-lengthweighted.cpp -o exe/astral-lengthweighted_precise.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-hybrid.cpp -o exe/astral-hybrid.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-hybrid.cpp -o exe/astral-hybrid_precise.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/master-site.cpp -o exe/master-site.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/master-pair.cpp -o exe/master-pair.exe