# Makefile for dopsearch

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn -I./lib -I/data/fit/aj/scan -I/rst/code/src.lib/analysis/aacgm/aacgm.1.10/src/ -I/rst/code/src.lib/analysis/ribmath/ -I$(IPATH)/analysis -I/usr/include
OBJS=dopsearch.o dopsearch_lib.o
SRC=dopsearch.c dopsearch_lib.c
OUTPUT = ./dopsearch
LIBS=-loldfit.1 -lgtabw.1 -loldgtabw.1 -lrpos.1 -lcfit.1 -lfit.1 -lrscan.1 -lradar.1 -ldmap.1 -lopt.1 -lscan2.1 -lrtime.1 -lrcnv.1 -laacgm.1 -ligrf.1 -lribmath.1
SLIB=-lm
include $(MAKEBIN)
