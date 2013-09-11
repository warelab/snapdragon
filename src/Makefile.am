if MACOS
    suff = -mt
endif
LDADD = -lboost_serialization$(suff)
if BUILD_FASTBIT
   EXTRA_BUILD=fastbit
endif
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++

SUBDIRS = bvec preprocessing kmerizer util programs $(EXTRA_BUILD)
lib_LTLIBRARIES = libsnapdragon.la
libsnapdragon_la_SOURCES =
nodist_EXTRA_libsnapdragon_la_SOURCES = dummy.cxx
libsnapdragon_la_LIBADD = \
	util/libutil.la \
	bvec/libbvec.la \
	kmerizer/libkmerizer.la
