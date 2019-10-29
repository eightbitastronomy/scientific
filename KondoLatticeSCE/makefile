#not in current use: -DUSEMAGCASE32 -DCOMMENSURATE_Q -DCRUDELUMP
#         debugging: -DDEBUG -DDEBUG1 -DDEBUG2 -DDEBUG3
#           current: -DUSECASEJ32  -DKSUM -DLUMP -DDSFMT_MEXP=521 -DSINGLEPROC
CC=cc -DLUMP -DKSUM -DDSFMT_MEXP=521
LFLAGS=-lm -lgsl -lgslcblas
DEPS=kdefs.h bare.h fio.h nint.h sm.h kzerocore.h kfinitecore.h sdwfinitecore.h sdwzerocore.h krandom.h bulkcore.h quartic.h


all: sdwzero sdwfinite sdwfinite-snc kondozero kondofinite sdwcorrect-m sdwcorrect-s kondocorrect-s kondoflatfinite
####################################### CURRENT CODE #################################################################################
afmbulk: afm_bulk.o bulkcore.o nint.o sm.o fio.o
	$(CC) -o afmbulk afm_bulk.o bulkcore.o nint.o sm.o fio.o $(LFLAGS)
kondobulk: kondo_bulk.o bulkcore.o nint.o sm.o fio.o krandom.o dSFMT/dSFMT.o
	$(CC) -o kondobulk kondo_bulk.o bulkcore.o nint.o sm.o fio.o krandom.o dSFMT/dSFMT.o $(LFLAGS)
afmisozero: afmiso_zero.o nint.o sm.o fio.o bare.o afmisozerocore.o quartic.o
	$(CC) -o afmisozero afmiso_zero.o nint.o sm.o fio.o bare.o afmisozerocore.o quartic.o $(LFLAGS)
afmisofinite: afmiso_finite.o nint.o sm.o fio.o afmisofinitecore.o
	$(CC) -o afmisofinite afmiso_finite.o nint.o sm.o fio.o afmisofinitecore.o $(LFLAGS)
kondozero: kondo_zero_general.o fio.o nint.o kzerocore.o sm.o bare.o
	$(CC) -o kondozero kondo_zero_general.o fio.o nint.o kzerocore.o sm.o bare.o $(LFLAGS)
kondofinite: kondo_finite_general.o fio.o nint.o kfinitecore.o sm.o bare.o
	$(CC) -o kondofinite kondo_finite_general.o fio.o nint.o kfinitecore.o sm.o bare.o $(LFLAGS)
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $<
clean:
	rm *.o

