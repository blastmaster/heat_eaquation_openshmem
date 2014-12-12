CC=oshcc
SRC_ROOT=src/

SOURCE= \
		$(SRC_ROOT)shmem_heat2d.c \
		$(SRC_ROOT)sharedblk.c

OUT=heat

CCFLAGS=-DOSH

all:
	$(CC) $(CCFLAGS) $(SOURCE) -o $(OUT)

clean:
	rm $(OUT)
