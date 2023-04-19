CXX = clang++

all:
	$(CXX) main.cc vec3.cc utils.cc -I./include -o rt -O2 -std=c++20 -Wall -Werror

debug:
	$(CXX) main.cc vec3.cc utils.cc -I./include -o rt -O0 -g -std=c++20 -Wall -Werror
