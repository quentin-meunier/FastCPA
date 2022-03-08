

CFLAGS = -O2

all: fastcpa cpa_ref

fastcpa: fastcpa.o common.o
	gcc -o $@ $^ -lm

cpa_ref: cpa_ref.o common.o
	gcc -o $@ $^ -lm

common.o: common.c common.h
	gcc -std=c99 -c $< -o $@ $(CFLAGS)

fastcpa.o: fastcpa.c common.h
	gcc -std=c99 -c $< -o $@ $(CFLAGS)

cpa_ref.o: cpa_ref.c common.h
	gcc -std=c99 -c $< -o $@ $(CFLAGS)



