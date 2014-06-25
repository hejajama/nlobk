CXXFLAGS = `gsl-config --cflags` -g -pedantic -I ../amplitudelib  
LDFLAGS = `gsl-config --libs` -lm 

include filelist.m

all: rbk

rbk: $(OBJECTS) 
	g++ $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) ../amplitudelib/libamplitude.a -o nlobk 
.cpp.o:
	 g++ $(CXXFLAGS) $< -c -o $@
.c.o:
	gcc $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJECTS) $(AMPLITUDELIBO)
	rm -f nlobk	
