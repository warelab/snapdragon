check_PROGRAMS = test-kmerizer test-bvec test-freqmap
TESTS = histo.test $(check_PROGRAMS)
AM_DEFAULT_SOURCE_EXT = .cpp
CC = g++
SOURCES = test-kmerizer.cpp test-bvec.cpp test-freqmap.cpp test.h
ACLOCAL_AMFLAGS = -I m4
GTEST_ROOT = $(top_srcdir)/vendor/googletest
AM_CXXFLAGS = -I$(top_srcdir)/src -I$(GTEST_ROOT)/include
AM_LDFLAGS = -L$(top_srcdir)/src -L$(GTEST_ROOT)/lib
if MACOS
    suff = -mt
endif
LDADD = -lsnapdragon \
	-lboost_serialization$(suff) \
	-lboost_thread$(suff) \
	-lboost_system$(suff) \
	-lgtest
EXTRA_DIST = $(TESTS)