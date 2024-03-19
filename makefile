all: dir astral astral-pro astral-hybrid caster-site caster-site_branchlength caster-pair waster-site
	echo "*** Installation complete! ***"

mac: dir astral astral-pro astral-hybrid
	echo "*** Installation complete! ***"

astral: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral.cpp -o bin/astral || g++ -std=c++17 -O2 -pthread src/astral.cpp -o bin/astral
	g++ -std=gnu++11 -march=native -D CASTLES -Ofast -pthread src/astral.cpp -o bin/astral4_experimental

astral_int128: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral.cpp -o bin/astral_int128
	
astral-pro: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro || g++ -std=c++17 -O2 -pthread src/astral-pro.cpp -o bin/astral-pro
	g++ -std=gnu++11 -march=native -D CASTLES -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro3_experimental
	
astral-pro_int128: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro_int128

astral-weighted: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-weighted.cpp -o bin/astral-weighted
	
astral-weighted_precise: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-weighted.cpp -o bin/astral-weighted_precise
	
astral-lengthweighted: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-lengthweighted.cpp -o bin/astral-lengthweighted
	
astral-lengthweighted_precise: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-lengthweighted.cpp -o bin/astral-lengthweighted_precise
	
astral-hybrid: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-hybrid.cpp -o bin/astral-hybrid || g++ -std=c++17 -O2 -pthread src/astral-hybrid.cpp -o bin/astral-hybrid
	
astral-hybrid_precise: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-hybrid.cpp -o bin/astral-hybrid_precise
	
caster-site: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/caster-site.cpp -o bin/caster-site || g++ -std=gnu++17 -O2 -pthread src/caster-site.cpp -o bin/caster-site

caster-site_branchlength: dir
	g++ -std=gnu++11 -march=native -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -Ofast -pthread src/caster-site.cpp -o bin/caster-site_branchlength || g++ -std=gnu++17 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread src/caster-site.cpp -o bin/caster-site_branchlength_experimental
	
caster-pair: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/caster-pair.cpp -o bin/caster-pair || g++ -std=gnu++17 -O2 -pthread src/caster-pair.cpp -o bin/caster-pair

waster-site: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/waster-site.cpp -o bin/waster-site || g++ -std=gnu++17 -O2 -pthread src/waster-site.cpp -o bin/waster-site

sister: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/sister.cpp -o bin/sister || g++ -std=gnu++17 -O2 -pthread src/sister.cpp -o bin/sister
	
dir:
	g++ -v 2>&1 | tail -n 1
	echo 'If installation failed, please ensure g++ version >= 7.5.0'
	mkdir -p bin
	
clean:
	rm bin/*

aaa-gpu:
	nvcc -O3 -rdc=true -c src/biallelic-cuda.cu -o bin/biallelic-cuda.o
	nvcc -dlink -o bin/biallelic-cuda_link.o bin/biallelic-cuda.o -lcudadevrt -lcudart
	g++ -std=gnu++11 -march=native -Ofast -pthread -D"USE_CUDA" bin/biallelic-cuda.o bin/biallelic-cuda_link.o src/asterisk-biallelic.cpp -o bin/asterisk-biallelic-cuda -L/usr/local/cuda/lib64 -lcudart -lcudadevrt
	rm bin/*.o
	
tutorial: all
	bin/astral -H
	bin/astral-pro -H
	bin/astral-hybrid -H
	bin/caster-site -H
	bin/waster-site -H
