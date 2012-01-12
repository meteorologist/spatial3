# 
# $Id: configure.ac Tue Oct 18 2011 aleksandarb $
#

AC_PREREQ(2.59)
AC_INIT([spatial3], [0.1.0], [aleksandarb@met.no])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AM_INIT_AUTOMAKE([no-dependencies])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_LN_S
#AC_PROG_INSTALL

AC_SUBST([SPATIAL3_SO_VERSION],  [0:0:0])
AC_SUBST([SPATIAL3_API_VERSION], [0.1.0])


METNO_PROG_PKG_CONFIG

#-----------------------------------------------------------------------------
# Checks
#-----------------------------------------------------------------------------

#
# C++
#
AC_LANG_PUSH([C++])

AC_CXX_HAVE_STD

AC_LANG_POP([C++])

#
#METNO_REQUIRE_HDF5
#METNO_REQUIRE_BOOST_FEATURE([date_time])
METNO_REQUIRE_GEOS
METNO_REQUIRE_PROJ
#METNO_REQUIRE_MYSQL
#METNO_REQUIRE_MYSQLPP
METNO_REQUIRE_BOOST_MT


# Output

METNO_PROGRAM_SUFFIX

AC_CONFIG_FILES([
    Makefile
    src/Makefile
	tests/Makefile
])

AC_OUTPUT
