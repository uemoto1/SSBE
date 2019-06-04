include Makefile.inc

TARGET = SSBE 

OBJS = \
misc/unusedvar.o \
parallel/salmon_communication.o \
parallel/salmon_parallel.o \
common/structures.o \
common/pack_unpack.o \
main.o


.SUFFIXES: .f90

%.o: %.f90
	$(FC) -c -o $@ $^ $(FLAGS)

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FLAGS) $(LIBS)

.PHONY: all clean 

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJS) *.mod


