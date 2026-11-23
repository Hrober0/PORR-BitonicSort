# Sekwencyjne sortowanie bitoniczne.

CC      := gcc
CFLAGS  := -std=c11 -O2 -Wall -Wextra -Icommon
LDLIBS  := -lm

COMMON_SRC := common/bitonic_common.c
COMMON_HDR := common/bitonic_common.h

.PHONY: all seq clean calibrate

all: seq

seq: bitonic_seq

bitonic_seq: seq/bitonic_seq.c $(COMMON_SRC) $(COMMON_HDR)
	$(CC) $(CFLAGS) -o $@ seq/bitonic_seq.c $(COMMON_SRC) $(LDLIBS)

clean:
	rm -f bitonic_seq

# Kalibracja log2_n od 20 do 28. Cel: około 120 s.
calibrate: bitonic_seq
	@mkdir -p results
	@echo "Sequential calibration: log2_n from 20 to 28 (target ~120 s)."
	@echo "Rows are appended to results/calibrate.csv (run id = log2_n)."
	@for k in 20 21 22 23 24 25 26 27 28; do \
		echo "---- log2_n=$$k ----"; \
		./bitonic_seq $$k 1 1 results/calibrate.csv $$k || exit $$?; \
	done
	@echo "Choose n_max as the log2_n whose time_s is closest to 120 seconds."
