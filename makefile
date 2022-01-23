all: astral astral-pro astral-weighted astral-lengthweighted astral-hybrid asterisk

astral:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral.cpp -o bin/astral

astral-pro: 
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro

astral-weighted:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-weighted.cpp -o bin/astral-weighted
	
astral-lengthweighted:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-lengthweighted.cpp -o bin/astral-lengthweighted
	
astral-hybrid:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-hybrid.cpp -o bin/astral-hybrid
	
asterisk:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/asterisk.cpp -o bin/asterisk
	
clean:
	rm bin/*
