prefix ?= $(HOME)/opt
CC ?= gcc

PKGCONFIG := PKG_CONFIG_PATH="$(prefix)/lib/pkgconfig" pkg-config
COMPILE := $(CC) -std=c99 -Wall -Wextra -pedantic -fPIC -O3 -fopenmp -pipe -g -MMD -I../include `$(PKGCONFIG) --cflags mandelbrot-numerics cairo` -c
LINK    := $(CC) -fopenmp -shared -g `$(PKGCONFIG) --libs mandelbrot-numerics cairo`
LIBRARY := libmandelbrot-graphics
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))

all: $(LIBRARY).a $(LIBRARY).so pkgconfig/mandelbrot-graphics.pc

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(LIBRARY).a $(LIBRARY).so pkgconfig/mandelbrot-graphics.pc

install: $(LIBRARY).a $(LIBRARY).so ../include/mandelbrot-graphics.h pkgconfig/mandelbrot-graphics.pc 
	install -d "$(prefix)/include" "$(prefix)/lib" "$(prefix)/lib/pkgconfig"
	install -m 644 -t "$(prefix)/include" ../include/mandelbrot-graphics.h
	install -m 644 -t "$(prefix)/lib" $(LIBRARY).a $(LIBRARY).so
	install -m 644 -t "$(prefix)/lib/pkgconfig" pkgconfig/mandelbrot-graphics.pc

$(LIBRARY).a: $(OBJECTS)
	@echo "A       $@" ; ar -rs $@ $^ || ( echo "ERROR   ar -rs $@ $^" && false )

$(LIBRARY).so: $(OBJECTS)
	@echo "SO      $@" ; $(LINK) -o $@ $^ -lmpc -lmpfr -lgmp -lm || ( echo "ERROR   $(LINK) -o $@ $^ -lmpc -lmpfr -lgmp -lm" && false )

%.o: %.c
	@echo "O       $@" ; $(COMPILE) -o $@ $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

pkgconfig/mandelbrot-graphics.pc: pkgconfig/mandelbrot-graphics.pc.in
	@echo "PC      $@" ; ( echo "prefix=$(prefix)" ; cat pkgconfig/mandelbrot-graphics.pc.in ) > pkgconfig/mandelbrot-graphics.pc || ( echo 'ERROR    ( echo "prefix=$(prefix)" ; cat pkgconfig/mandelbrot-graphics.pc.in ) > pkgconfig/mandelbrot-graphics.pc' && false )

.SUFFIXES:
.PHONY: all clean install

-include $(DEPENDS)
