prefix ?= $(HOME)/opt
CC ?= gcc

PKGCONFIG := PKG_CONFIG_PATH="$(prefix)/lib/pkgconfig" pkg-config
COMPILE := $(CC) -g -std=c99 -Wall -Wextra -fPIC -O3 -pipe -MMD `$(PKGCONFIG) --cflags mandelbrot-graphics mandelbrot-numerics mandelbrot-symbolics`
LIBS    := `$(PKGCONFIG) --libs mandelbrot-graphics mandelbrot-numerics mandelbrot-symbolics` -lm -lgmp
OBJECTS := $(patsubst %.c,%.o,$(wildcard *.c))
DEPENDS := $(patsubst %.o,%.d,$(OBJECTS))
EXES    := $(patsubst %.o,%,$(OBJECTS))

all: $(EXES)

clean:
	@echo "CLEAN" ; rm -f $(OBJECTS) $(DEPENDS) $(EXES)

install: $(EXES)
	install -d "$(prefix)/bin"
	install -m 755 -t "$(prefix)/bin" $(EXES)

%: %.o
	@echo "EXE     $@" ; gcc -o $@ $< $(LIBS) || ( echo "ERROR   gcc -o $@ $< $(LIBS)" && false )

%.o: %.c
	@echo "O       $@" ; $(COMPILE) -o $@ -c $< || ( echo "ERROR   $(COMPILE) -o $@ $<" && false )

.SUFFIXES:
.PHONY: all clean install
.SECONDARY: $(OBJECTS)

-include $(DEPENDS)
