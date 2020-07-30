EXTRA_CXXFLAGS=
CXXFLAGS=-O3 -Wall -std=c++17 $(EXTRA_CXXFLAGS)

main: pjpeg.o pjpeg_compress.o pjpeg_decompress.o
	g++ -o pjpeg pjpeg.o pjpeg_compress.o pjpeg_decompress.o

pjpeg.o: pjpeg.cpp
	g++ $(CXXFLAGS) -c pjpeg.cpp

pjpeg_compress.o: pjpeg_compress.cpp
	g++ $(CXXFLAGS) -c pjpeg_compress.cpp

pjpeg_decompress.o: pjpeg_decompress.cpp
	g++ $(CXXFLAGS) -c pjpeg_decompress.cpp

clean:
	rm -f pjpeg *.o
