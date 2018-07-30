CXXFLAGS = `gsl-config --cflags` -O2 -fopenmp -pedantic -I ../amplitudelib_v2/ -I ../amplitudelib_v2/tools/ -I src/
LDFLAGS = `gsl-config --libs` -lm 

include filelist.m

all: rbk

rbk: $(OBJECTS) src/main.o 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) src/main.o ../amplitudelib_v2/build/lib/libamplitude.a -o nlobk 
tools: tools/nlobklibtest.o $(OBJECTS)
	${CXX} $(CXXFLAGS) -I src/ $(LDFLAGS) $(OBJECTS)  tools/nlobklibtest.o ../amplitudelib_v2/build/lib/libamplitude.a -o tools/nlobklibtest
.cpp.o: 
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO)
	rm -f nlobk	
