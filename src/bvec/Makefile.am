noinst_LTLIBRARIES = libbvec.la
libbvec_la_SOURCES = bvec.cpp bvec32.cpp bvec.h bvec32.h kseq.h
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
if MACOS
   suff = -mt
endif
LDADD = -lboost_serialization$(suff)
