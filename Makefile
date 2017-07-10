CXXFLAGS = `gsl-config --cflags` -O2 -fopenmp -pedantic -I ../amplitudelib_v2/amplitudelib2/  
LDFLAGS = `gsl-config --libs` -lm 

include filelist.m

all: rbk

rbk: $(OBJECTS) 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) ../amplitudelib_v2/amplitudelib2/libamplitude.a -o nlobk 
tools: tools/nlobklibtest.o $(OBJECTS)
	${CXX} $(CXXFLAGS) -I src/ $(LDFLAGS) $(OBJECTS)  tools/nlobklibtest.o ../amplitudelib_v2/libamplitude.a -o tools/nlobklibtest
.cpp.o: 
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@
clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO)
	rm -f nlobk	
