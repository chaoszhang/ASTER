mkdir exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral.cpp -o exe/astral.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral.cpp -o exe/astral_int128.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-pro.cpp -o exe/astral-pro.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-pro.cpp -o exe/astral-pro_int128.exe
::g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-weighted.cpp -o exe/astral-weighted.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-weighted.cpp -o exe/astral-weighted_precise.exe
::g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-lengthweighted.cpp -o exe/astral-lengthweighted.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-lengthweighted.cpp -o exe/astral-lengthweighted_precise.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-hybrid.cpp -o exe/astral-hybrid.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-hybrid.cpp -o exe/astral-hybrid_precise.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/caster-site.cpp -o exe/caster-site.exe
g++ -std=gnu++11 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread -static -static-libgcc -static-libstdc++ src/caster-site.cpp -o exe/caster-site_branchlength.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/caster-pair.cpp -o exe/caster-pair.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/waster-site.cpp -o exe/waster-site.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/sister.cpp -o exe/sister.exe