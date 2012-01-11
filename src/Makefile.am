## using Automake and Libtool to build a shared library

spatialite_INCLUDEPATHS = -I$(srcdir)/backend -I$(srcdir)/backend/headers

geographiclib_INCLUDEPATHS = -I$(srcdir)/backend/GeographicLib

spatialite_CPPFLAGS = -D_LARGE_FILE=1 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE=1 -DOMIT_FREEXL=1 \
                      -DSQLITE_ENABLE_RTREE=1 -DNDEBUG=1 -DSPL_AMALGAMATION=1 -DGEOS_ADVANCED=1 \
                      -Wall -Wextra
#            @OMIT_GEOS_FLAGS@ @GEOSADVANCED_FLAGS@ @OMIT_PROJ_FLAGS@ \
#            @OMIT_ICONV_FLAGS@ @OMIT_MATHSQL_FLAGS@ \
#            @OMIT_EPSG_FLAGS@ @OMIT_GEOCALLBACKS_FLAGS@ @OMIT_FREEXL_FLAGS@

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
#### AUTOMAKE_OPTIONS = subdir-objects

## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
#### ACLOCAL_AMFLAGS = -I m4

## Define a libtool archive target "libmetgm-@SPATIAL3_API_VERSION@.la", with
## @SPATIAL3_API_VERSION@ substituted into the generated Makefile at configure
## time.
## The libtool archive file (.la) will be installed into the directory named
## by the predefined variable $(bindir), along with the actual shared library
## file (.so).
lib_LTLIBRARIES = libspatial3.la

## Define the source file list for the "libmetgm-@metgm_API_VERSION@.la"
## target.  Note that @SPATIAL3_API_VERSION@ is not interpreted by Automake and
## will therefore be treated as if it were literally part of the target name,
## and the variable name derived from that.
## The file extension .cc is recognized by Automake, and makes it produce
## rules which invoke the C++ compiler to produce a libtool object file (.lo)
## from each source file.  Note that it is not necessary to list header files
## which are already listed elsewhere in a _HEADERS variable assignment.
libspatial3_la_SOURCES = \
                       backend/headers/spatialite/gaiaexif.h \
                       backend/headers/spatialite/gaiaaux.h \
                       backend/headers/spatialite/gaiageo.h \
                       backend/headers/spatialite/gg_const.h \
                       backend/headers/spatialite/gg_structs.h \
                       backend/headers/spatialite/gg_core.h \
                       backend/headers/spatialite/gg_mbr.h \
                       backend/headers/spatialite/gg_formats.h \
                       backend/headers/spatialite/gg_dynamic.h \
                       backend/headers/spatialite/gg_advanced.h \
                       backend/headers/spatialite/sqlite3.h \
                       backend/headers/spatialite/sqlite3ext.h \
                       backend/headers/spatialite/spatialite.h \
                       backend/headers/spatialite.h \
                       backend/GeographicLib/Accumulator.hpp \
                       backend/GeographicLib/AlbersEqualArea.hpp \
                       backend/GeographicLib/AzimuthalEquidistant.hpp \
                       backend/GeographicLib/CassiniSoldner.hpp \
                       backend/GeographicLib/CircularEngine.hpp \
                       backend/GeographicLib/Constants.hpp \
                       backend/GeographicLib/DMS.hpp \
                       backend/GeographicLib/EllipticFunction.hpp \
                       backend/GeographicLib/GeoCoords.hpp \
                       backend/GeographicLib/Geocentric.hpp \
                       backend/GeographicLib/Geodesic.hpp \
                       backend/GeographicLib/GeodesicLine.hpp \
                       backend/GeographicLib/Geoid.hpp \
                       backend/GeographicLib/Gnomonic.hpp \
                       backend/GeographicLib/GravityCircle.hpp \
                       backend/GeographicLib/GravityModel.hpp \
                       backend/GeographicLib/LambertConformalConic.hpp \
                       backend/GeographicLib/LocalCartesian.hpp \
                       backend/GeographicLib/MGRS.hpp \
                       backend/GeographicLib/MagneticCircle.hpp \
                       backend/GeographicLib/MagneticModel.hpp \
                       backend/GeographicLib/Math.hpp \
                       backend/GeographicLib/NormalGravity.hpp \
                       backend/GeographicLib/OSGB.hpp \
                       backend/GeographicLib/PolarStereographic.hpp \
                       backend/GeographicLib/PolygonArea.hpp \
                       backend/GeographicLib/SphericalEngine.hpp \
                       backend/GeographicLib/SphericalHarmonic.hpp \
                       backend/GeographicLib/SphericalHarmonic1.hpp \
                       backend/GeographicLib/SphericalHarmonic2.hpp \
                       backend/GeographicLib/TransverseMercator.hpp \
                       backend/GeographicLib/TransverseMercatorExact.hpp \
                       backend/GeographicLib/UTMUPS.hpp \
                       backend/GeographicLib/Utility.hpp \
                       backend/GeographicLib/Config.h \
                       backend/AlbersEqualArea.cpp \
                       backend/AzimuthalEquidistant.cpp \
                       backend/CassiniSoldner.cpp \
                       backend/CircularEngine.cpp \
                       backend/DMS.cpp \
                       backend/EllipticFunction.cpp \
                       backend/GeoCoords.cpp \
                       backend/Geocentric.cpp \
                       backend/Geodesic.cpp \
                       backend/GeodesicLine.cpp \
                       backend/Geoid.cpp \
                       backend/Gnomonic.cpp \
                       backend/GravityCircle.cpp \
                       backend/GravityModel.cpp \
                       backend/LambertConformalConic.cpp \
                       backend/LocalCartesian.cpp \
                       backend/MGRS.cpp \
                       backend/MagneticCircle.cpp \
                       backend/MagneticModel.cpp \
                       backend/NormalGravity.cpp \
                       backend/OSGB.cpp \
                       backend/PolarStereographic.cpp \
                       backend/PolygonArea.cpp \
                       backend/SphericalEngine.cpp \
                       backend/TransverseMercator.cpp \
                       backend/TransverseMercatorExact.cpp \
                       backend/UTMUPS.cpp \
                       backend/Utility.cpp \
                       backend/spatialite.c \
                       backend/sqlite3.c \
					   geometryimpl.hpp \
                       coordinate.hpp \
					   coordinate.cpp \
                       geometry.hpp \
                       geometry.cpp \
                       geometryd2.hpp \
                       geometryd2.cpp \
                       geometryll.hpp \
                       geometryll.cpp \
                       region.hpp

libspatial3_la_CPPFLAGS = \
                        ${spatialite_INCLUDEPATHS} \
                        ${geographiclib_INCLUDEPATHS} \
                        ${GEOS_CPPFLAGS} \
                        ${PROJ_CPPFLAGS} \
                        ${spatialite_CPPFLAGS}
                        


## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
libspatial3_la_LDFLAGS = \
                       -version-info $(SPATIAL3_SO_VERSION) \
                       ${GEOS_LDFLAGS} ${GEOS_LIBS} \
                       ${PROJ_LDFLAGS} ${PROJ_LIBS}

## Define the list of public header files and their install location.  The
## nobase_ prefix instructs Automake to not strip the directory part from each
## filename, in order to avoid the need to define separate file lists for each
## installation directory.  This only works if the directory hierarchy in the
## source tree matches the hierarchy at the install location, however.
spatial3_includedir = $(includedir)
dist_spatial3_include_HEADERS = \
                              coordinate.hpp \
                              geometry.hpp   \
                              geometryd2.hpp \
                              geometryll.hpp \
                              region.hpp

## The generated configuration header is installed in its own subdirectory of
## $(libdir).  The reason for this is that the configuration information put
## into this header file describes the target platform the installed library
## has been built for.  Thus the file must not be installed into a location
## intended for architecture-independent files, as defined by the Filesystem
## Hierarchy Standard (FHS).
## The nodist_ prefix instructs Automake to not generate rules for including
## the listed files in the distribution on 'make dist'.  Files that are listed
## in _HEADERS variables are normally included in the distribution, but the
## configuration header file is generated at configure time and should not be
## shipped with the source tarball.
spatial3_libincludedir = $(libdir)/include

## Install the generated pkg-config file (.pc) into the expected location for
## architecture-dependent package configuration information.  Occasionally,
## pkg-config files are also used for architecture-independent data packages,
## in which case the correct install location would be $(datadir)/pkgconfig.
### pkgconfigdir = $(libdir)/pkgconfig
### pkgconfig_DATA = metgm.pc
