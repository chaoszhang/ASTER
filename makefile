all: astral astral-pro

astral:
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral.cpp -o bin/astral

astral-pro: 
	g++ -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro
	
clean:
	rm bin/*
