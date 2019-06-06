include Makefile.inc

TARGET = SSBE 

OBJS = \
misc/unusedvar.o \
io/salmon_file.o \
parallel/salmon_communication.o \
parallel/salmon_parallel.o \
common/structures.o \
common/pack_unpack.o \
input_parameter.o \
sbe.o \
main.o

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FLAGS) $(LIBS)

.SUFFIXES: .f90

%.o: %.f90
	$(FC) -c -o $@ $^ $(FLAGS)

.PHONY: all clean 

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJS) *.mod

input_parameter.f90: input_parameter.py
	python input_parameter.py > $@
