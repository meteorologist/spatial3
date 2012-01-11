#AM_CFLAGS = -D_LARGE_FILE=1 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE=1 \
#            -DSQLITE_ENABLE_RTREE=1 -DNDEBUG=1 -DSPL_AMALGAMATION=1 \
#	        @OMIT_GEOS_FLAGS@ @GEOSADVANCED_FLAGS@ @OMIT_PROJ_FLAGS@ \
#            @OMIT_ICONV_FLAGS@ @OMIT_MATHSQL_FLAGS@ \
#	        @OMIT_EPSG_FLAGS@ @OMIT_GEOCALLBACKS_FLAGS@ @OMIT_FREEXL_FLAGS@

#INCLUDES = @CFLAGS@
#INCLUDES += -I$(top_srcdir)/src/backend/headers

#SUBDIRS = headers

sources = \
        headers/spatialite.h \
        headers/spatialite/gaiaexif.h \
	headers/spatialite/gaiaaux.h \
	headers/spatialite/gaiageo.h \
	headers/spatialite/gg_const.h \
	headers/spatialite/gg_structs.h \
	headers/spatialite/gg_core.h \
	headers/spatialite/gg_mbr.h \
	headers/spatialite/gg_formats.h \
	headers/spatialite/gg_dynamic.h \
	headers/spatialite/gg_advanced.h \
	headers/spatialite/sqlite3.h \
	headers/spatialite/sqlite3ext.h \
	headers/spatialite/spatialite \
        spatialite.c \
	sqlite3.c
