mkdir exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral.cpp -o exe/astral4.exe
copy exe\astral4.exe exe\astral.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral.cpp -o exe/astral_int128.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-pro.cpp -o exe/astral-pro3.exe
copy exe\astral-pro3.exe exe\astral-pro.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-pro.cpp -o exe/astral-pro_int128.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-hybrid.cpp -o exe/wastral.exe
copy exe\wastral.exe exe\astral-hybrid.exe
::g++ -std=gnu++11 -D LARGE_DATA -O2 -pthread -static -static-libgcc -static-libstdc++ src/astral-hybrid.cpp -o exe/astral-hybrid_precise.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/caster-site.cpp -o exe/caster-site.exe
g++ -std=gnu++11 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread -static -static-libgcc -static-libstdc++ src/caster-site.cpp -o exe/caster-site_branchlength_experimental.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/caster-pair.cpp -o exe/caster-pair.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/waster-site.cpp -o exe/waster-site.exe
g++ -std=gnu++11 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread -static -static-libgcc -static-libstdc++ src/waster-site.cpp -o exe/waster-site_branchlength_experimental.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/waster-ng.cpp -o exe/waster-ng_experimental.exe
g++ -std=gnu++11 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread -static -static-libgcc -static-libstdc++ src/waster-ng.cpp -o exe/waster-ng_branchlength_experimental.exe
g++ -std=gnu++11 -O2 -pthread -static -static-libgcc -static-libstdc++ src/sister.cpp -o exe/sister.exe