all : docks
docks: docks.cpp
	g++ -O3 -std=c++11 -Wall -o docks docks.cpp
