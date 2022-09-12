all: dir astral astral_int128 astral-pro astral-pro_int128 astral-weighted astral-weighted_precise astral-lengthweighted astral-lengthweighted_precise astral-hybrid astral-hybrid_precise asterisk asterisk-hky
	echo "*** Installation complete! ***"

astral:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral.cpp -o bin/astral

astral_int128:
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral.cpp -o bin/astral_int128
	
astral-pro: 
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro
	
astral-pro_int128:
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro_int128

astral-weighted:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-weighted.cpp -o bin/astral-weighted
	
astral-weighted_precise:
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-weighted.cpp -o bin/astral-weighted_precise
	
astral-lengthweighted:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-lengthweighted.cpp -o bin/astral-lengthweighted
	
astral-lengthweighted_precise:
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-lengthweighted.cpp -o bin/astral-lengthweighted_precise
	
astral-hybrid:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-hybrid.cpp -o bin/astral-hybrid
	
astral-hybrid_precise:
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-hybrid.cpp -o bin/astral-hybrid_precise
	
asterisk:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/asterisk.cpp -o bin/asterisk

asterisk-hky:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/asterisk-hky.cpp -o bin/asterisk-hky
	
dir:
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
