all: dir astral astral_coalescent_unit astral-pro wastral caster-site caster-site_branchlength caster-pair waster waster_branchlength
	echo "*** Installation complete! ***"

mac: dir astral astral-pro wastral
	echo "*** Installation complete! ***"

astral: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral.cpp -o bin/astral4 || g++ -std=c++17 -O2 -pthread src/astral.cpp -o bin/astral4
	cp bin/astral4 bin/astral

astral_int128: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral.cpp -o bin/astral4_int128

astral_coalescent_unit:	
	g++ -std=gnu++11 -march=native -D COALESCENT_UNIT -Ofast -pthread src/astral.cpp -o bin/astral4_coalescent_unit || g++ -std=c++17 -D COALESCENT_UNIT -O2 -pthread src/astral.cpp -o bin/astral4_coalescent_unit

astral-pro: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro3 || g++ -std=c++17 -O2 -pthread src/astral-pro.cpp -o bin/astral-pro3
	cp bin/astral-pro3 bin/astral-pro
	
astral-pro_int128: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro3_int128

astral-pro_coalescent_unit: dir
	g++ -std=gnu++11 -march=native -D COALESCENT_UNIT -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro3_coalescent_unit || g++ -std=c++17 -D COALESCENT_UNIT -O2 -pthread src/astral-pro.cpp -o bin/astral-pro3_coalescent_unit

wastral: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-hybrid.cpp -o bin/wastral || g++ -std=c++17 -O2 -pthread src/astral-hybrid.cpp -o bin/wastral
	cp bin/wastral bin/astral-hybrid
	
wastral_precise: dir
	g++ -std=gnu++11 -march=native -D LARGE_DATA -Ofast -pthread src/astral-hybrid.cpp -o bin/wastral_precise
	
caster-site: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/caster-site.cpp -o bin/caster-site || g++ -std=gnu++17 -O2 -pthread src/caster-site.cpp -o bin/caster-site

caster-site_branchlength: dir
	g++ -std=gnu++11 -march=native -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -Ofast -pthread src/caster-site.cpp -o bin/caster-site_branchlength || g++ -std=gnu++17 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread src/caster-site.cpp -o bin/caster-site_branchlength_experimental
	
caster-pair: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/caster-pair.cpp -o bin/caster-pair || g++ -std=gnu++17 -O2 -pthread src/caster-pair.cpp -o bin/caster-pair

waster-old_branchlength: dir
	g++ -std=gnu++11 -march=native -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -Ofast -pthread src/waster-old.cpp -o bin/waster-old_branchlength || g++ -std=gnu++17 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread src/waster-old.cpp -o bin/waster-old_branchlength_experimental
	
waster-old: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/waster-old.cpp -o bin/waster-old || g++ -std=gnu++17 -O2 -pthread src/waster-old.cpp -o bin/waster-old

waster_branchlength: dir
	g++ -std=gnu++11 -march=native -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -Ofast -pthread src/waster.cpp -o bin/waster_branchlength || g++ -std=gnu++17 -D CUSTOMIZED_ANNOTATION_TERMINAL_LENGTH -O2 -pthread src/waster.cpp -o bin/waster_branchlength
	
waster: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/waster.cpp -o bin/waster || g++ -std=gnu++17 -O2 -pthread src/waster.cpp -o bin/waster

sister: dir
	g++ -std=gnu++11 -march=native -Ofast -pthread src/sister.cpp -o bin/sister || g++ -std=gnu++17 -O2 -pthread src/sister.cpp -o bin/sister
	
dir:
	g++ -v 2>&1 | tail -n 1
	echo 'If installation failed, please ensure g++ version >= 7.5.0'
	mkdir -p bin
	
clean:
	rm bin/*
	
tutorial: all
	mkdir -p tutorial
	bin/astral4 -H
	bin/astral-pro3 -H
	bin/wastral -H
	bin/caster-site -H
	bin/waster-site -H
