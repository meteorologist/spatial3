bin_PROGRAMS = simpletest

simpletest_SOURCES = simpletest.cpp

spatial3_ldflags = -L$(PWD)/../src/.libs
spatial3_libs = -lspatial3
libs = -lpthread -ldl

spatial3_cppflags= -I$(PWD)/../src

simpletest_CPPFLAGS = \
                     ${GEOS_CPPFLAGS} \
                     ${PROJ_CPPFLAGS} \
                     ${spatial3_cppflags}

simpletest_LDFLAGS = \
                    ${libs} \
                    ${GEOS_LDFLAGS} ${GEOS_LIBS} \
                    ${PROJ_LDFLAGS} ${PROJ_LIBS} \
					${spatial3_ldflags} ${spatial3_libs}
