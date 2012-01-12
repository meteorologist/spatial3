
#ifndef S3_GEOMETRY_IMPL_H_
#define S3_GEOMETRY_IMPL_H_

#include "coordinate.hpp"

// GEOC C++
//
#include <geos/geom/Geometry.h>
#include <geos/geom/Point.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/operation/IsSimpleOp.h>
#include <geos/operation/buffer/BufferParameters.h>
#include <geos/operation/valid/IsValidOp.h>
#include <geos/operation/union/UnaryUnionOp.h>
#include <geos/operation/polygonize/Polygonizer.h>
#include <geos/simplify/DouglasPeuckerSimplifier.h>
#include <geos/simplify/TopologyPreservingSimplifier.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>
#include <geos/io/ParseException.h>
#include <geos/util/UniqueCoordinateArrayFilter.h>

// boost
//
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

namespace gg = geos::geom;
namespace gi = geos::io;
namespace go = geos::operation;

namespace metno { namespace s3 {

// ## -------------------------------------------------------------------------------------------------------
//     static GEOS related utilities
// ## -------------------------------------------------------------------------------------------------------
    static void destroy(gg::Geometry* g)
    {
        if(!g) return;

        const gg::GeometryFactory* factory = g->getFactory();
        factory->destroyGeometry(g);
    }

    // caller must take ownership of returned geometry
    static gg::Geometry* add_point_linestring(gg::Geometry const* g, gg::Coordinate c)
    {
        gg::CoordinateSequence* coords = g->getCoordinates();
        coords->add(c);
        return g->getFactory()->createLineString(coords);
    }

    // caller must take ownership of returned geometry
    static gg::Geometry* close2d_ring(gg::Geometry const* g)
    {
        gg::LineString const* ring = dynamic_cast<gg::LineString const*>(g);
        gg::Geometry* closed;

        // close the ring if not already closed
        if(!ring->isClosed()) {
            gg::Coordinate first(ring->getCoordinate()->x, ring->getCoordinate()->y);
            closed = add_point_linestring(ring, first);
        } else {
            closed = ring->clone();
        }
        return closed;
    }

    static gg::Geometry* make_geos_friendly_ring(gg::Geometry const* r)
    {
        gg::Geometry* closed = close2d_ring(r);

        // return 0 for collapsed ring (after closeup)
        while(closed->getNumPoints() < 4) {
            gg::Geometry* tmp = add_point_linestring(closed, gg::Coordinate(closed->getCoordinate()->x, closed->getCoordinate()->y));
            tmp->getFactory()->destroyGeometry(closed);
            closed = tmp;
        }
        return closed;
    }

    static gg::Geometry* make_geos_friendly_polygon(gg::Geometry* p)
    {
        gg::Polygon* poly = dynamic_cast<gg::Polygon*>(p);

        // If the polygon has no rings there's nothing to do
        if(poly->getNumInteriorRing() == 0)
            return poly->clone();

        std::vector<gg::Geometry*>* holes = new std::vector<gg::Geometry*>();
        for(size_t i = 0; i < poly->getNumInteriorRing(); ++i) {
            gg::Geometry* closed = make_geos_friendly_ring(poly->getInteriorRingN(i));
            holes->push_back(closed);
        }

        assert(poly->getExteriorRing()->isClosed());

        gg::LinearRing* shell = p->getFactory()->createLinearRing(poly->getExteriorRing()->getCoordinates());

        return p->getFactory()->createPolygon(shell, holes);
    }

    /*
     * Fully node given linework
     */
    static gg::Geometry* node_lines(gg::Geometry const* lines)
    {
        gg::LineString const* ls = dynamic_cast<gg::LineString const*>(lines);
        gg::Geometry* noded;
        gg::Geometry* point;

        /*
         * Union with first geometry point, obtaining full noding
         * and dissolving of duplicated repeated points
         *
         * TODO: substitute this with UnaryUnion?
         */
        point = ls->getPointN(0);
        if(!point)
            return 0;

        noded = ls->Union(point);

        ls->getFactory()->destroyGeometry(point);

        if(!noded) {
            return 0;
        }

        return noded;
    }

    static gg::Geometry* extract_unique_points(gg::Geometry const* g)
    {
        namespace gu = geos::util;

        try
        {
            /* 1: extract points */
            std::vector<gg::Coordinate const*> coords;
            gu::UniqueCoordinateArrayFilter filter(coords);
            g->apply_ro(&filter);

            /* 2: for each point, create a geometry and put into a vector */
            std::vector<gg::Geometry*>* points = new std::vector<gg::Geometry*>();
            points->reserve(coords.size());
            const gg::GeometryFactory* factory = g->getFactory();
            std::vector<gg::Coordinate const*>::iterator it, itE;
            for (it = coords.begin(); it != coords.end(); ++it)
            {
                gg::Geometry* point = factory->createPoint(*(*it));
                points->push_back(point);
            }

            /* 3: create a multipoint */
            return g->getFactory()->createMultiPoint(points);

        }
        catch (const std::exception &e)
        {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "<<e.what()<<std::endl;
            return 0;
        }
        catch (...)
        {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : unknown exception throw "<<std::endl;;
            return 0;
        }
    }

    static gg::Geometry* polygonize(const gg::Geometry* const* vgeoms, size_t ngeoms)
    {
        gg::Geometry *out = 0;

        if(!vgeoms || !vgeoms[0])
            return out;

//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" vgeoms[0] WKT "<<vgeoms[0]->toString()<<std::endl;

        try
        {
            // Polygonize
            geos::operation::polygonize::Polygonizer polygonizer;
            for (size_t i = 0; i < ngeoms; ++i)
            {
                polygonizer.add(vgeoms[i]);
                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" polygonizer adding "<<vgeoms[0]->getGeometryType()<<std::endl;
            }

            std::vector<gg::Polygon*> *polys = polygonizer.getPolygons();

            assert(0 != polys);

            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" polygonizer returned #polygons "<<polys->size()<<std::endl;

            std::vector<gg::Geometry*> *polyvec = new std::vector<gg::Geometry *>(polys->size());

            for (size_t i = 0; i < polys->size(); ++i)
            {
                (*polyvec)[i] = (*polys)[i];
            }
            delete polys;
            polys = 0;

            const gg::GeometryFactory *gf = vgeoms[0]->getFactory();

            // The below takes ownership of the passed vector,
            // so we must *not* delete it
            out = gf->createGeometryCollection(polyvec);

//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" out WKT "<<out->toString()<<std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "<<e.what()<<std::endl;
            return 0;
        }
        catch (...)
        {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : Unknown exception thrown "<<std::endl;
            return 0;
        }

        return out;
    }

    static gg::Geometry* build_area(gg::Geometry const* geom_in)
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" geom_in TYPE "<<geom_in->getGeometryType()<<std::endl;
        gg::Geometry* shp;
        gg::Geometry const* vgeoms[1];
        size_t i;

        vgeoms[0] = geom_in;

        gg::Geometry* geos_result = polygonize(vgeoms, 1);
        if (!geos_result)
            return 0;

        /*
         * We should now have a collection
         */
        if(geos_result->getGeometryTypeId() != gg::GEOS_GEOMETRYCOLLECTION) {
            destroy(geos_result);
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" Unexpected return from GEOSpolygonize "<<std::endl;
            return 0;
        }

        size_t ngeoms = geos_result->getNumGeometries();

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" geos_result TYPE "<<geos_result->getGeometryType()<<std::endl;

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" geos_result->getNumGeometries() " << ngeoms <<std::endl;

        /*
         * No geometries in collection, early out
         */
        if(ngeoms == 0) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" EARLY EXIT "<<std::endl;
            return geos_result;
        }

        /*
         * Return first geometry if we only have one in collection,
         * to avoid the unnecessary Geometry clone below.
         */
        if(ngeoms == 1)
        {
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            const gg::Geometry* tmp = geos_result->getGeometryN(0);

            if(!tmp) {
                geom_in->getFactory()->destroyGeometry(geos_result);
                return 0;
            }
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            shp = tmp->clone();
            geom_in->getFactory()->destroyGeometry(geos_result); /* only safe after the clone above */
            return shp;
        }

        /*
         * Iteratively invoke symdifference on outer rings
         * as suggested by Carl Anderson:
         * postgis-devel/2005-December/001805.html
         */
        shp = 0;
        for (i=0; i < ngeoms; ++i)
        {
            /*
             * Construct a Polygon from geometry i exterior ring
             * We don't use GEOSGeom_clone on the ExteriorRing
             * due to a bug in CAPI contained in GEOS 2.2 branch
             * failing to properly return a LinearRing from
             * a LinearRing clone.
             */
            gg::Polygon const* geomN = dynamic_cast<gg::Polygon const*>(geos_result->getGeometryN(i));
            gg::LineString const* extls = geomN->getExteriorRing();
            gg::CoordinateSequence* sq = extls->getCoordinates();

            gg::Geometry* extring =
                    geom_in->getFactory()->createPolygon(geom_in->getFactory()->createLinearRing(sq), 0);

            if (extring == 0) /* exception */
            {
//                lwerror("GEOSCreatePolygon threw an exception");
                return 0;
            }

            if(shp == 0) {
                shp = extring;
//                LWDEBUGF(3, "GEOSpolygonize: shp:%s",
//                         lwgeom_to_ewkt(GEOS2LWGEOM(shp, 0)));
            } else {
                gg::Geometry* tmp = shp->symDifference(extring);
//                LWDEBUGF(3, "GEOSpolygonize: SymDifference(%s, %s):%s",
//                         lwgeom_to_ewkt(GEOS2LWGEOM(shp, 0)),
//                         lwgeom_to_ewkt(GEOS2LWGEOM(extring, 0)),
//                         lwgeom_to_ewkt(GEOS2LWGEOM(tmp, 0))
//                         );
                geom_in->getFactory()->destroyGeometry(shp);
                geom_in->getFactory()->destroyGeometry(extring);
                shp = tmp;
            }
        }

        geom_in->getFactory()->destroyGeometry(geos_result);

        return shp;
    }

    static gg::Geometry* make_valid_polygon(gg::Geometry const* p)
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" p TYPE is "<<p->getGeometryType() <<std::endl;
        gg::Geometry* gout;
        gg::Geometry* collapse_points;
        std::vector<gg::Geometry*> vgeoms; /* One for area, one for cut-edges */
        size_t nvgeoms = 0;

        assert(p->getGeometryTypeId() == gg::GEOS_POLYGON || p->getGeometryTypeId() == gg::GEOS_MULTIPOLYGON);

        gg::Geometry* bound = p->getBoundary();
        if(!bound)
            return 0;

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" bound TYPE is "<<bound->getGeometryType()<<std::endl;

        gg::Geometry* cut_edges = node_lines(bound);
        if(!cut_edges) {
            p->getFactory()->destroyGeometry(bound);
            return 0;
        }

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" cut_edges TYPE is "<<cut_edges->getGeometryType()<<std::endl;
        /*
         * NOTE: the noding process may drop lines collapsing to points.
         * We want to retrive any of those
         */
        {
            gg::Geometry* pi;
            gg::Geometry* po;

            pi = extract_unique_points(bound);
            if (pi == 0) {
                p->getFactory()->destroyGeometry(bound);
                return 0;
            }

//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" bound unique points "<<pi->toString()<<std::endl;

            po = extract_unique_points(cut_edges);
            if (po == 0) {
                p->getFactory()->destroyGeometry(bound);
                p->getFactory()->destroyGeometry(pi);
                return 0;
            }
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" cut_edges unique points "<<po->toString()<<std::endl;

            collapse_points = pi->difference(po);
            if(collapse_points == 0) {
                p->getFactory()->destroyGeometry(bound);
                p->getFactory()->destroyGeometry(pi);
                p->getFactory()->destroyGeometry(po);
                return 0;
            }
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" collapse unique points "<<collapse_points->toString()<<std::endl;


            p->getFactory()->destroyGeometry(pi);
            p->getFactory()->destroyGeometry(po);
        }
        p->getFactory()->destroyGeometry(bound);

        // And use an empty geometry as initial "area"
        gg::Geometry* area = p->getFactory()->createPolygon();
        if(!area) {
            p->getFactory()->destroyGeometry(cut_edges);
            return 0;
        }

        /*
         * See if an area can be build with the remaining edges
         * and if it can, symdifference with the original area.
         * Iterate this until no more polygons can be created
         * with left-over edges.
         */
        while(cut_edges->getNumGeometries())
        {
            gg::Geometry* new_area=0;
            gg::Geometry* new_area_bound=0;
            gg::Geometry* symdif=0;
            gg::Geometry* new_cut_edges=0;

            /*
             * ASSUMPTION: cut_edges should already be fully noded
             */

            new_area = build_area(cut_edges);
            if (!new_area)   /* must be an exception */
            {
                p->getFactory()->destroyGeometry(cut_edges);
                p->getFactory()->destroyGeometry(area);
//                lwnotice("LWGEOM_GEOS_buildArea() threw an error: %s", lwgeom_geos_errmsg);
                return 0;
            }


            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" new_area "<<new_area->toString()<<std::endl;

            if (new_area->isEmpty()){
                /* no more rings can be build with thes edges */
                p->getFactory()->destroyGeometry(new_area);
                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" BREAK "<<std::endl;
                break;
            }

            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" We succeeded in building a ring ! "<<std::endl;
            /*
             * Save the new ring boundaries first (to compute
             * further cut edges later)
             */
            new_area_bound = new_area->getBoundary();
            if(!new_area_bound)
            {
                /*
                 * We did check for empty area already so
                 * this must be some other error
                */
                p->getFactory()->destroyGeometry(new_area);
                p->getFactory()->destroyGeometry(area);
                return 0;
            }

            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" new area "<<new_area->toString()<<std::endl;
            /*
             * Now symdif new and old area
             */
            symdif = area->symDifference(new_area);

//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" SYM "<<symdif->toString()<<std::endl;

            if (!symdif)   /* must be an exception */
            {
                p->getFactory()->destroyGeometry(cut_edges);
                p->getFactory()->destroyGeometry(new_area);
                p->getFactory()->destroyGeometry(new_area_bound);
                p->getFactory()->destroyGeometry(area);
//                lwnotice("GEOSSymDifference() threw an error: %s", lwgeom_geos_errmsg);
                return 0;
            }

            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;

            p->getFactory()->destroyGeometry(area);
            p->getFactory()->destroyGeometry(new_area);
            area = symdif;
            symdif = 0;
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" AREA "<<area->toString()<<std::endl;
            /*
             * Now let's re-set geos_cut_edges with what's left
             * from the original boundary.
             * ASSUMPTION: only the previous cut-edges can be
             *             left, so we don't need to reconsider
             *             the whole original boundaries
             */
            new_cut_edges = cut_edges->difference(new_area_bound);
            p->getFactory()->destroyGeometry(new_area_bound);
            if ( ! new_cut_edges )   /* an exception ? */
            {
                /* cleanup and throw */
                p->getFactory()->destroyGeometry(cut_edges);
                p->getFactory()->destroyGeometry(area);
//                lwnotice("GEOSDifference() threw an error: %s", lwgeom_geos_errmsg);
                return NULL;
            }
            p->getFactory()->destroyGeometry(cut_edges);
            cut_edges = new_cut_edges;
        }

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" area isEmpty "<<(area->isEmpty()?"TRUE":"FALSE") <<std::endl;

        if(!area->isEmpty()) {
            vgeoms.push_back(area);
            nvgeoms++;
        } else {
            p->getFactory()->destroyGeometry(area);
        }

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" cut_edges isEmpty "<<(cut_edges->isEmpty()?"TRUE":"FALSE") <<std::endl;

        if(!cut_edges->isEmpty()) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            vgeoms.push_back(cut_edges);
            nvgeoms++;
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            p->getFactory()->destroyGeometry(cut_edges);
        }

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" collapse_points isEmpty "<<(collapse_points->isEmpty()?"TRUE":"FALSE") <<std::endl;

        if(!collapse_points->isEmpty()) {
            vgeoms.push_back(collapse_points);
            nvgeoms++;
        } else {
            p->getFactory()->destroyGeometry(collapse_points);
        }
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;

        if (nvgeoms == 1) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" nvgeoms == 1 " <<std::endl;
            /* Return cut edges */
            gout = vgeoms.at(0);
        } else {
            /* Collect areas and lines (if any line) */
            gout = p->getFactory()->createGeometryCollection(vgeoms);
            if (!gout)   /* an exception again */
            {
                /* cleanup and throw */
//                lwnotice("GEOSGeom_createCollection() threw an error: %s", lwgeom_geos_errmsg);
                /* TODO: cleanup! */
                return 0;
            }
        }

//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" gout "<<gout->toString()<<std::endl;
        return gout;
    }

    static gg::Geometry* make_clean_polygon(gg::Geometry const* poly_in)
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;

        gg::Geometry* gout;

        gout = make_valid_polygon(poly_in);

        if(!gout) {
            return 0;
        }

        /* Check dimensionality is the same as input */
        if(poly_in->getDimension() != gout->getDimension()) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__
                     <<" dimensional collapse "
                     <<poly_in->getDimension()
                     <<" to "
                     <<gout->getDimension()
                     <<std::endl;
            return 0;
        }

        /* Check that the output is not a collection if the input wasn't */
        if(gout->getGeometryTypeId() == gg::GEOS_GEOMETRYCOLLECTION
             && poly_in->getGeometryTypeId() != gg::GEOS_GEOMETRYCOLLECTION)
        {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__
                     <<" mixed-type output "
                     <<gout->getGeometryType()
                     <<" from single-type input "
                     <<poly_in->getGeometryType()
                     <<std::endl;
            return 0;
        }

        // TODO: -- maybe --
        /* gout := ST_ForceRHR(gout); */
        /* gout = ST_RemoveRepeatedPoints(gout); */

        return gout;
    }

    static gg::Geometry* line_split_by_line(const gg::Geometry* line_in, const gg::Geometry* blade_in)
    {
        gg::Geometry* gdiff; /* difference */
        gg::Geometry* g1 = line_in->clone();
        gg::Geometry* g2 = blade_in->clone();

        /* Possible outcomes:
         *
         *  1. The lines do not cross or overlap
         *      -> Return a collection with single element
         *  2. The lines cross
         *      -> Return a collection of all elements resulting from the split
         */

        /* If interior intersecton is linear we can't split */
        if(g1->relate(g2, std::string("1********"))) {
            destroy(g1);
            destroy(g2);
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" Splitter line has linear intersection with input "<<std::endl;
            return 0;
        }

        gdiff = g1->difference(g2);
        destroy(g1);
        destroy(g2);
        if(gdiff == 0) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" empty difference with splitter "<<std::endl;
            return 0;
        }

        if(gdiff->getGeometryTypeId() == gg::GEOS_GEOMETRYCOLLECTION) {
            // collection already
            return gdiff;
        } else {
            // have to make collection
            std::vector<gg::Geometry*>* pGeoms(new std::vector<gg::Geometry*>());
            pGeoms->push_back(gdiff);
            return gdiff->getFactory()->createGeometryCollection(pGeoms);
        }
    }

    static gg::Geometry* polygon_split_by_line(const gg::Geometry* polygon, const gg::Geometry* line)
    {
        gg::Geometry* polygon_ = polygon->clone();
        gg::Geometry* line_ = line->clone();

        /* Possible outcomes:
         *
         *  1. The line does not split the polygon
         *      -> Return a collection with single element
         *  2. The line does split the polygon
         *      -> Return a collection of all elements resulting from the split
         */

        gg::Geometry* polygon_bounds_ = polygon_->getBoundary();
        if(!polygon_bounds_) {
            destroy(polygon_);
            return 0;
        }

        gg::Geometry* vgeoms = polygon_bounds_->Union(line_);
        if(!vgeoms) {
            destroy(polygon_);
            destroy(polygon_bounds_);
            destroy(line_);
            return 0;
        }

        geos::operation::polygonize::Polygonizer polygonizer;
        polygonizer.add(vgeoms);
        std::vector<gg::Polygon*>*  polys = polygonizer.getPolygons();

        if(polys->size() < 2) {
            destroy(polygon_);
            destroy(polygon_bounds_);
            destroy(line_);
            destroy(vgeoms);
            delete polys;
            return 0;
        }

        std::vector<gg::Geometry*>* polyvec(new std::vector<gg::Geometry*>(polys->size()));
        for(size_t i = 0; i < polys->size(); ++i) (*polyvec)[i] = (*polys)[i];
        delete polys;

        const gg::GeometryFactory *gf = polygon_->getFactory();

        gg::Geometry* polygons = gf->createGeometryCollection(polyvec);
        if(!polygons) {
            destroy(polygon_);
            destroy(polygon_bounds_);
            destroy(line_);
            destroy(vgeoms);
            return 0;
        }

        assert(polygons->getGeometryTypeId() == gg::GEOS_GEOMETRYCOLLECTION);

    #if PARANOIA_LEVEL > 0
            if ( GEOSGeometryTypeId(polygons) != COLLECTIONTYPE )
            {
                    GEOSGeom_destroy(g1);
                    GEOSGeom_destroy(g2);
                    GEOSGeom_destroy(g1_bounds);
                    GEOSGeom_destroy((GEOSGeometry*)vgeoms[0]);
                    GEOSGeom_destroy(polygons);
                    lwerror("Unexpected return from GEOSpolygonize");
                    return 0;
            }
    #endif

        /* We should now have all polygons, just skip
         * the ones which are in holes of the original
         * geometries and return the rest in a collection
         */
        size_t n = polygons->getNumGeometries();
        std::vector<gg::Geometry*>* geometries(new std::vector<gg::Geometry*>());
        for (size_t i = 0; i < n; ++i)
        {
            gg::Geometry* pos; /* point on surface */
            const gg::Geometry* p = polygons->getGeometryN(i);

            pos = p->getInteriorPoint();
            if(!pos) {
                destroy(polygon_);
                destroy(polygon_bounds_);
                destroy(line_);
                destroy(vgeoms);
                destroy(polygons);
                return 0;
            }

            if(polygon_->contains(pos))
                geometries->push_back(p->clone());

            destroy(pos);

        }

        destroy(polygon_);
        destroy(polygon_bounds_);
        destroy(line_);
        destroy(vgeoms);
        destroy(polygons);

        return gf->createGeometryCollection(geometries);
    }

    static gg::Geometry* collection_split_by_line(const gg::Geometry* coll_in, const gg::Geometry* blade_in)
    {
        std::vector<gg::Geometry*>* parts(new std::vector<gg::Geometry*>());

        const gg::GeometryFactory* gf = coll_in->getFactory();

        gg::Geometry* l = gf->createEmptyGeometry();
        gg::Geometry* r = gf->createEmptyGeometry();

        for(size_t i = 0; i < coll_in->getNumGeometries(); ++i) {
            gg::Geometry* split = polygon_split_by_line(coll_in->getGeometryN(i), blade_in);
            if(!split) {
                return 0;
            } else {
//                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" WKT "<<split->toString()<<std::endl;
                gg::Geometry* tmpL = l->Union(split->getGeometryN(0));
                gg::Geometry* tmpR = r->Union(split->getGeometryN(1));
                destroy(split);
                destroy(l);
                destroy(r);
                l = tmpL;
                r = tmpR;
            }
        }

        parts->push_back(l);
        parts->push_back(r);

        return coll_in->getFactory()->createGeometryCollection(parts);
    }

// ## -------------------------------------------------------------------------------------------------------
//     EnvelopeImpl
// ## -------------------------------------------------------------------------------------------------------
        class EnvelopeImpl {
        public:
            EnvelopeImpl()
                : envelope_(new gg::Envelope())
            {}
            EnvelopeImpl(double x1, double x2, double y1, double y2)
                : envelope_(new gg::Envelope(x1, x2, y1, y2))
            {}

            EnvelopeImpl(const Coordinate &p1, const Coordinate &p2)
                : envelope_(new gg::Envelope(gg::Coordinate(p1.getX(), p1.getY()), gg::Coordinate(p2.getX(), p2.getY())))
            {}

            double getArea() const {
                return envelope_->getArea();
            }

            double getMaxX() const {
                return envelope_->getMaxX();
            }

            double getMaxY() const {
                return envelope_->getMaxY();
            }
            double getMinX() const {
                return envelope_->getMinX();
            }

            double getMinY() const {
                return envelope_->getMinY();
            }

            std::string toString() const {
                return envelope_->toString();
            }

            bool centre(Coordinate &centre) const {
                gg::Coordinate ggC;
                bool result = envelope_->centre(ggC);
                centre.setX(ggC.x);
                centre.setY(ggC.y);
                return result;
            }

            void setToNull() {
                envelope_->setToNull();
            }

            bool isNull() const {
                return envelope_->isNull();
            }

            bool equals(const EnvelopeImpl *other) const {
                return envelope_->equals(other->envelope_.get());
            }

            double getWidth() const {
                return envelope_->getWidth();
            }

            double getHeight() const {
                return envelope_->getHeight();
            }

            bool intersection(const EnvelopeImpl &env, EnvelopeImpl &result) const {
                return envelope_->intersection(*env.envelope_, *result.envelope_);
            }

            bool intersects(EnvelopeImpl const* other) const {
                return envelope_->intersects(other->envelope_.get());
            }

            bool intersects(Coordinate const& c) const {
                return envelope_->intersects(gg::Coordinate(c.getX(), c.getY()));
            }

        protected:
        private:
            boost::shared_ptr<gg::Envelope> envelope_;
        };

// ## -------------------------------------------------------------------------------------------------------
//     static GEOS related utilities
// ## -------------------------------------------------------------------------------------------------------
    template <int SRID>
    class GeometryImpl {
    public:

        GeometryImpl()
        {
            pGeometry_ = boost::shared_ptr<gg::Geometry>(pFactory_->createEmptyGeometry(), destroy);
        }

        GeometryImpl(gg::Geometry* geom)
        {
            pGeometry_ = boost::shared_ptr<gg::Geometry>(geom, destroy);
        }

//        GeometryImpl(GeometryImpl const &rhs)
//        {
////            pGeometry_ = rhs.pGeometry_;
//        }

//        GeometryImpl(GeometryImpl const* g)
//        {
////            pGeometry_ = g->pGeometry_;
//        }

//        GeometryImpl &operator=(GeometryImpl const &rhs)
//        {
//             if (this == &rhs)
//                 return *this;

////            pGeometry_ = rhs.pGeometry_;
//        }

        virtual ~GeometryImpl() { }

        virtual std::string toWKT() const;

        virtual int          getSRID() const;
        virtual GeometryType getGeometryType() const;

        virtual bool isEmpty()     const;
        virtual bool isValid()     const;
        virtual bool isRectangle() const;
        virtual bool isSimple()    const;

        virtual CoordinatePtr getNonSimpleLocation() const;

//        virtual double getArea() const;
//        virtual double getLength() const;
        virtual size_t getNumVertices() const;
        virtual size_t getNumGeometries() const;
        virtual void getCoordinates(CoordinatePtrSequence& coordinates) const;
        virtual Coordinate getCentroid() const;
        EnvelopeImpl getEnvelope() const;
        virtual Coordinate getLowerLeftCorner() const;
        virtual Coordinate getUpperRightCorner() const;
        virtual boost::shared_ptr<GeometryImpl<SRID> > getGeometryN(size_t n) const;
        DimensionType getDimension() const;

        // releation operators
        virtual bool intersects(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool touches(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool crosses(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool contains(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool within(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool equals(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool overlaps(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool covers(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual bool coveredBy(boost::shared_ptr<GeometryImpl> const& g) const;

        // operators
        virtual boost::shared_ptr<GeometryImpl> Intersection(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual boost::shared_ptr<GeometryImpl> Difference(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual boost::shared_ptr<GeometryImpl> Union(boost::shared_ptr<GeometryImpl> const& g) const;
        virtual boost::shared_ptr<GeometryImpl<SRID> > Buffer(double distance) const;

        // creation
        static boost::shared_ptr<GeometryImpl> createPoint(double x = 0.0, double y = 0.0);
        static boost::shared_ptr<GeometryImpl> createPolygon(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryImpl> createLineString(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryImpl> createLinearRing(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryImpl> createMultiPoint(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryImpl> createMultiLineString(std::vector<CoordinatePtrSequence> const& coordinates);
        static boost::shared_ptr<GeometryImpl> createMultiPolygon(std::vector<CoordinatePtrSequence> const& coordinates);
        static boost::shared_ptr<GeometryImpl> fromWKT(const std::string& wkt);

        // simplification
        virtual void simplifyDP(double tolerance);
        virtual void simplifyTP(double tolerance);

        // advanced functionality
        boost::shared_ptr<GeometryImpl> makeGeosFriendly() const;
        virtual boost::shared_ptr<GeometryImpl> makeValid() const;
        virtual boost::shared_ptr<GeometryImpl> cleanGeometry() const;
        boost::shared_ptr<GeometryImpl> splitBy(boost::shared_ptr<GeometryImpl> const& g) const;

    protected:

        static boost::shared_ptr<GeometryImpl> fromRaw(gg::Geometry* gg);

        static boost::shared_ptr<gg::PrecisionModel>  pPM_;
        static boost::shared_ptr<gg::GeometryFactory> pFactory_;
        static boost::shared_ptr<gi::WKTReader>       pWtkReader_;
        static boost::shared_ptr<gi::WKTWriter>       pWtkWriter_;


        boost::shared_ptr<gg::Geometry>               pGeometry_;
    };

    template<int SRID>
    class GeometryCollectionImpl : public GeometryImpl<SRID> {
    public:
        GeometryCollectionImpl() {
             this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createGeometryCollection(), destroy);
        }
    };

    template<int SRID>
    class PointImpl : virtual public GeometryImpl<SRID> {
    public:

        PointImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createPoint(), destroy);
        }
        ~PointImpl() { }

        double getX();
        double getY();
        void   setX(double x);
        void   setY(double y);

    private:

    };

    template<int SRID>
    class LineStringImpl : public GeometryImpl<SRID> {
    public:
        LineStringImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createLineString());
        }

        LineStringImpl(gg::Geometry* geom) : GeometryImpl<SRID>(geom) { }

        bool isRing()   const;
        bool isClosed() const;

        boost::shared_ptr<PointImpl<SRID> > getEndPoint() const;
        boost::shared_ptr<PointImpl<SRID> > getStartPoint() const;

//        boost::shared_ptr<Geometry> reverse();

    protected:

    };

    template<int SRID>
    class LinearRingImpl : public LineStringImpl<SRID> {
    public:
        LinearRingImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createLinearRing());
        }
    };

    template<int SRID>
    class PolygonImpl : public GeometryImpl<SRID> {
    public:
        PolygonImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createPolygon());
        }

        boost::shared_ptr<GeometryImpl<SRID> > getExteriorRing();

        // advanced functionality
        boost::shared_ptr<GeometryImpl<SRID> > makeValid() const;
        boost::shared_ptr<GeometryImpl<SRID> > cleanGeometry() const;
    protected:

    };

    template <int SRID>
    class MultiPolygonImpl : public GeometryCollectionImpl<SRID> {
    public:
        MultiPolygonImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createMultiPolygon());
        }

        boost::shared_ptr<GeometryImpl<SRID> > makeValid() const;
        boost::shared_ptr<GeometryImpl<SRID> > cleanGeometry() const;
    };

    template <int SRID>
    class MultiPointImpl : public GeometryCollectionImpl<SRID> {
    public:
        MultiPointImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createMultiPoint());
        }
    };

    template <int SRID>
    class MultiLineStringImpl : public GeometryCollectionImpl<SRID> {
    public:
        MultiLineStringImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createMultiLineString());
        }
    };


    template <int SRID>
    boost::shared_ptr<gg::PrecisionModel> GeometryImpl<SRID>::pPM_(new gg::PrecisionModel());

    template <int SRID>
    boost::shared_ptr<gg::GeometryFactory> GeometryImpl<SRID>::pFactory_(new gg::GeometryFactory(pPM_.get(), SRID));

    template <int SRID>
    boost::shared_ptr<gi::WKTReader> GeometryImpl<SRID>::pWtkReader_(new gi::WKTReader(pFactory_.get()));

    template <int SRID>
    boost::shared_ptr<gi::WKTWriter> GeometryImpl<SRID>::pWtkWriter_(new gi::WKTWriter());

    template <int SRID>
    std::string GeometryImpl<SRID>::toWKT() const
    {
        if(pGeometry_.get())
            return pGeometry_->toString();
        else
            return std::string();
    }

    template <int SRID>
    int GeometryImpl<SRID>::getSRID() const
    {
        return GeometryImpl<SRID>::pFactory_->getSRID();
    }

    template <int SRID>
    GeometryType GeometryImpl<SRID>::getGeometryType() const
    {
        if(!this->pGeometry_.get()) {
            return S3_EMPTY;
        } else {
            switch(this->pGeometry_->getGeometryTypeId()) {
            case gg::GEOS_POINT:
                return S3_POINT;
            case gg::GEOS_LINESTRING:
                return S3_LINESTRING;
            case gg::GEOS_LINEARRING:
                return S3_LINEARRING;
            case gg::GEOS_POLYGON:
                return S3_POLYGON;
            case gg::GEOS_MULTIPOINT:
                return S3_MULTIPOINT;
            case gg::GEOS_MULTILINESTRING:
                return S3_MULTILINESTRING;
            case gg::GEOS_MULTIPOLYGON:
                return S3_MULTIPOLYGON;
            case gg::GEOS_GEOMETRYCOLLECTION:
                return S3_GEOMETRYCOLLECTION;
            }
        }
    }

    template <int SRID>
    bool GeometryImpl<SRID>::isEmpty() const
    {
        if(pGeometry_.get())
            return pGeometry_->isEmpty();
        else
            return true;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::isValid() const
    {
        if(pGeometry_.get()) {
            geos::operation::valid::IsValidOp validop(this->pGeometry_.get());
            if(!validop.isValid()){
                geos::operation::valid::TopologyValidationError* tve = validop.getValidationError();
                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" NON VALID GEOMETRY FOUND : "<<tve->toString()<<std::endl;
            }
            return pGeometry_->isValid();
        } else {
            return false;
        }
    }

    template <int SRID>
    bool GeometryImpl<SRID>::isSimple() const
    {
        if(pGeometry_.get()) {
            go::IsSimpleOp test(*pGeometry_.get());
            return test.isSimple();
        } else {
//            CoordinatePtr nonSimple = getNonSimpleLocation();
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" non-simple location "<< nonSimple->toString() <<std::endl;
            return false;
        }
    }

    template <int SRID>
    CoordinatePtr GeometryImpl<SRID>::getNonSimpleLocation() const
    {
        if(pGeometry_.get()) {
            go::IsSimpleOp test(*pGeometry_.get());

            // isSimple must be called first
            test.isSimple();

            const gg::Coordinate* location = test.getNonSimpleLocation();

            if(location) {
                return CoordinatePtr(new Coordinate(location->x, location->y));
            } else {
                return CoordinatePtr();
            }
        } else {
            return CoordinatePtr();
        }
    }

    template <int SRID>
    bool GeometryImpl<SRID>::isRectangle() const
    {
        if(pGeometry_.get())
            return pGeometry_->isRectangle();
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::equals(boost::shared_ptr<GeometryImpl> const& g) const
    {
        return this->pGeometry_->equals(g->pGeometry_.get());
    }

    template <int SRID>
    size_t GeometryImpl<SRID>::getNumVertices() const
    {
        if(this->pGeometry_.get())
            return this->pGeometry_->getNumPoints();
        else
            return 0;
    }

    template <int SRID>
    size_t GeometryImpl<SRID>::getNumGeometries() const
    {
        if(this->pGeometry_.get())
            return this->pGeometry_->getNumGeometries();
        else
            return 0;
    }

    template <int SRID>
    void GeometryImpl<SRID>::getCoordinates(CoordinatePtrSequence& coordinates) const
    {
        coordinates.clear();
        if(pGeometry_.get()) {
            gg::CoordinateSequence* sequence = pGeometry_->getCoordinates();
            if(sequence) {
                for(size_t i = 0; i < sequence->getSize(); ++i) {
                    coordinates.push_back(CoordinatePtr(new Coordinate(sequence->getAt(i).x, sequence->getAt(i).y)));
                }
            }
        }
    }

    template <int SRID>
    Coordinate GeometryImpl<SRID>::getCentroid() const
    {
        if(pGeometry_.get()) {
            gg::Coordinate ggC;
            pGeometry_->getCentroid(ggC);
            return Coordinate(ggC.x, ggC.y);
        } else {
            return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
    }

    template <int SRID>
    EnvelopeImpl GeometryImpl<SRID>::getEnvelope() const
    {
        if(pGeometry_.get()) {
            const gg::Envelope* ggEnv = pGeometry_->getEnvelopeInternal();
            return EnvelopeImpl(ggEnv->getMinX(), ggEnv->getMaxX(), ggEnv->getMinY(), ggEnv->getMaxY());
        } else {
            return EnvelopeImpl();
        }
    }

    template <int SRID>
    Coordinate GeometryImpl<SRID>::getUpperRightCorner() const
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
        if(pGeometry_.get()) {
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            if(pGeometry_->getEnvelopeInternal() != 0) {
//                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" ENVELOPE: "<<pGeometry_->getEnvelopeInternal()->toString()<<std::endl;
                return Coordinate(pGeometry_->getEnvelopeInternal()->getMaxX(), pGeometry_->getEnvelopeInternal()->getMaxY());
            }
        }

//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
        return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }

    template <int SRID>
    Coordinate GeometryImpl<SRID>::getLowerLeftCorner() const
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
        if(pGeometry_.get()) {
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            if(pGeometry_->getEnvelopeInternal() != 0) {
//                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" ENVELOPE: "<<pGeometry_->getEnvelopeInternal()->toString()<<std::endl;
                return Coordinate(pGeometry_->getEnvelopeInternal()->getMinX(), pGeometry_->getEnvelopeInternal()->getMinY());
            }
        }

//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
        return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::getGeometryN(size_t n) const
    {
        if(!this->pGeometry_.get() || pGeometry_->getGeometryN(n) == 0) {
            return boost::shared_ptr<GeometryImpl<SRID> >(new GeometryImpl<SRID>());
        } else {
            gg::Geometry* raw = pGeometry_->getGeometryN(n)->clone();
            return fromRaw(raw);
        }
    }

    template <int SRID>
    DimensionType GeometryImpl<SRID>::getDimension() const
    {
        switch(this->pGeometry_->getDimension()) {
        case gg::Dimension::A:
            return S3_DIM_SURFACE;
        case gg::Dimension::L:
            return S3_DIM_LINE;
        case gg::Dimension::P:
        default:
            return S3_DIM_POINT;
        }
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::Intersection(boost::shared_ptr<GeometryImpl<SRID> > const& g) const
    {
        gg::Geometry* ggRaw(pGeometry_->intersection(g->pGeometry_.get()));
        return fromRaw(ggRaw);
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::Difference(boost::shared_ptr<GeometryImpl<SRID> > const& g) const
    {
        gg::Geometry* ggRaw(pGeometry_->difference(g->pGeometry_.get()));
        return fromRaw(ggRaw);
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::Union(boost::shared_ptr<GeometryImpl<SRID> > const& g) const
    {
        gg::Geometry* ggRaw(pGeometry_->Union(g->pGeometry_.get()));
        return fromRaw(ggRaw);
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::Buffer(double distance) const
    {
        gg::Geometry* ggRaw(pGeometry_->buffer(distance, 2, geos::operation::buffer::BufferParameters::CAP_FLAT));
        return fromRaw(ggRaw);
    }

    template <int SRID>
    bool GeometryImpl<SRID>::intersects(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->intersects(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::contains(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->contains(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::touches(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->touches(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::crosses(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->crosses(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::within(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->within(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::overlaps(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->touches(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::covers(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->covers(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    bool GeometryImpl<SRID>::coveredBy(boost::shared_ptr<GeometryImpl> const& g) const
    {
        if(pGeometry_.get())
            return pGeometry_->coveredBy(g->pGeometry_.get());
        else
            return false;
    }

    template <int SRID>
    void GeometryImpl<SRID>::simplifyDP(double tolerance)
    {
        if(pGeometry_.get()) {
            geos::simplify::DouglasPeuckerSimplifier dps(pGeometry_.get());
            dps.setDistanceTolerance(tolerance);
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" before simplification "<<pGeometry_->toText()<<std::endl;
            std::auto_ptr<gg::Geometry> simplified = dps.getResultGeometry();
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" after simplification "<<simplified->toText()<<std::endl;
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" the tolerance was "<<tolerance<<std::endl;
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(simplified->clone());
        } else {
            return;
        }
    }

    template <int SRID>
    void GeometryImpl<SRID>::simplifyTP(double tolerance)
    {
        if(pGeometry_.get()) {
            geos::simplify::TopologyPreservingSimplifier tps(pGeometry_.get());
            tps.setDistanceTolerance(tolerance);
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" before simplification "<<pGeometry_->toText()<<std::endl;
            std::auto_ptr<gg::Geometry> simplified = tps.getResultGeometry();
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" after simplification "<<simplified->toText()<<std::endl;
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" the tolerance was "<<tolerance<<std::endl;
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(simplified->clone());
        } else {
            return;
        }
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::makeGeosFriendly() const
    {
        if(getGeometryType() == S3_POLYGON || getGeometryType() == S3_MULTIPOLYGON) {
            gg::Geometry* friendly = make_geos_friendly_polygon(this->pGeometry_.get());
            return GeometryImpl<SRID>::fromRaw(friendly);
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" not implemented for "<<pGeometry_->getGeometryType()<<" geometry type "<<std::endl;
            boost::shared_ptr<GeometryImpl<SRID> > empty(new GeometryImpl<SRID>());
            return empty;
        }
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::makeValid() const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" not implemented for "<<pGeometry_->getGeometryType()<<" geometry type "<<std::endl;
        boost::shared_ptr<GeometryImpl<SRID> > empty(new GeometryImpl<SRID>());
        return empty;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::cleanGeometry() const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" not implemented for "<<pGeometry_->getGeometryType()<<" geometry type "<<std::endl;
        boost::shared_ptr<GeometryImpl<SRID> > empty(new GeometryImpl<SRID>());
        return empty;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::splitBy(boost::shared_ptr<GeometryImpl<SRID> > const& g) const
    {
        if(getGeometryType() == S3_POLYGON && g->getGeometryType() == S3_LINESTRING) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            gg::Geometry* poly = this->pGeometry_.get();
            gg::Geometry* line = g->pGeometry_.get();
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            gg::Geometry* parts = polygon_split_by_line(poly, line);
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            return GeometryImpl<SRID>::fromRaw(parts);
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
        } else if(getGeometryType() == S3_MULTIPOLYGON && g->getGeometryType() == S3_LINESTRING) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            gg::Geometry* mpoly = this->pGeometry_.get();
            gg::Geometry* line = g->pGeometry_.get();
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
            gg::Geometry* parts = collection_split_by_line(mpoly, line);
            return GeometryImpl<SRID>::fromRaw(parts);
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" not implemented between "<<pGeometry_->getGeometryType()<<" and "<<g->pGeometry_->getGeometryType()<<std::endl;
            boost::shared_ptr<GeometryImpl<SRID> > empty(new GeometryImpl<SRID>());
            return empty;
        }
    }


// ## -------------------------------------------------------------------------------------------------------
//     PointImpl
// ## -------------------------------------------------------------------------------------------------------

    template <int SRID>
    double PointImpl<SRID>::getX()
    {
        if(this->pGeometry_->isEmpty())
            return std::numeric_limits<double>::quiet_NaN();
        gg::Point* pt = dynamic_cast<gg::Point*>(this->pGeometry_.get());
        return pt->getX();
    }

    template <int SRID>
    double PointImpl<SRID>::getY()
    {
        if(this->pGeometry_->isEmpty())
            return std::numeric_limits<double>::quiet_NaN();
        gg::Point* pt = dynamic_cast<gg::Point*>(this->pGeometry_.get());
        return pt->getY();
    }

    template <int SRID>
    void PointImpl<SRID>::setX(double newx)
    {
        double oldy = getY();
        this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createPoint(gg::Coordinate(newx, oldy)));
    }

    template <int SRID>
    void PointImpl<SRID>::setY(double newy)
    {
        double oldx = getX();
        this->pGeometry_ = boost::shared_ptr<gg::Point>(this->pFactory_->createPoint(gg::Coordinate(oldx, newy)));
    }

// ## -------------------------------------------------------------------------------------------------------
//     LineStringImpl
// ## -------------------------------------------------------------------------------------------------------

    template <int SRID>
    bool LineStringImpl<SRID>::isRing() const
    {
        gg::LineString* ls = dynamic_cast<gg::LineString*>(this->pGeometry_.get());
        return ls->isRing();
    }

    template <int SRID>
    bool LineStringImpl<SRID>::isClosed() const
    {
        gg::LineString* ls = dynamic_cast<gg::LineString*>(this->pGeometry_.get());
        return ls->isClosed();
    }

    template <int SRID>
    boost::shared_ptr<PointImpl<SRID> > LineStringImpl<SRID>::getEndPoint() const
    {
        gg::LineString* ls = dynamic_cast<gg::LineString*>(this->pGeometry_.get());
        std::auto_ptr<gg::Point> end(ls->getEndPoint());
        return boost::dynamic_pointer_cast<PointImpl<SRID> >(GeometryImpl<SRID>::createPoint(end->getX(), end->getY()));
    }

    template <int SRID>
    boost::shared_ptr<PointImpl<SRID> > LineStringImpl<SRID>::getStartPoint() const
    {
        gg::LineString* ls = dynamic_cast<gg::LineString*>(this->pGeometry_.get());
        std::auto_ptr<gg::Point> start(ls->getStartPoint());
        return boost::dynamic_pointer_cast<PointImpl<SRID> >(GeometryImpl<SRID>::createPoint(start->getX(), start->getY()));
    }


// ## -------------------------------------------------------------------------------------------------------
//     PolygonImpl
// ## -------------------------------------------------------------------------------------------------------

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > PolygonImpl<SRID>::getExteriorRing()
    {
        boost::shared_ptr<GeometryImpl<SRID> > empty(new LineStringImpl<SRID>());
        if(!this->pGeometry_.get() || this->pGeometry_->isEmpty())
            return empty;

        gg::Polygon* ggPolygon = dynamic_cast<gg::Polygon*>(this->pGeometry_.get());
        if(!ggPolygon)
            return empty;

        const gg::LineString* ggShell = ggPolygon->getExteriorRing();
        gg::CoordinateSequence* seq(ggShell->getCoordinates());
        gg::Geometry* gg = this->pFactory_->createLineString(seq);
        boost::shared_ptr<GeometryImpl<SRID> > ls(new LineStringImpl<SRID>(gg));
        return ls;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > PolygonImpl<SRID>::makeValid() const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
        gg::Geometry* friendly = make_geos_friendly_polygon(this->pGeometry_.get());
        gg::Geometry* valid = make_valid_polygon(friendly);
//        gg::Geometry* valid = make_valid_polygon_by_buffering_offending_point(friendly);

        destroy(friendly);

        return GeometryImpl<SRID>::fromRaw(valid);
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > PolygonImpl<SRID>::cleanGeometry() const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
        gg::Geometry* friendly = make_geos_friendly_polygon(this->pGeometry_.get());
        gg::Geometry* valid = make_clean_polygon(friendly);
        destroy(friendly);

        return GeometryImpl<SRID>::fromRaw(valid);
    }

    // ## -------------------------------------------------------------------------------------------------------
    //     MultiPolygonImpl
    // ## -------------------------------------------------------------------------------------------------------

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > MultiPolygonImpl<SRID>::makeValid() const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
        gg::Geometry* friendly = make_geos_friendly_polygon(this->pGeometry_.get());
        gg::Geometry* valid = make_valid_polygon(friendly);

        destroy(friendly);

        return GeometryImpl<SRID>::fromRaw(valid);
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > MultiPolygonImpl<SRID>::cleanGeometry() const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK " <<std::endl;
        gg::Geometry* friendly = make_geos_friendly_polygon(this->pGeometry_.get());
        gg::Geometry* valid = make_clean_polygon(friendly);
        destroy(friendly);

        return GeometryImpl<SRID>::fromRaw(valid);
    }

// ## -------------------------------------------------------------------------------------------------------
//     creation static functions
// ## -------------------------------------------------------------------------------------------------------

    template<int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::fromRaw(gg::Geometry* ggRaw)
    {
        if(!ggRaw) {
            boost::shared_ptr<GeometryImpl<SRID> > empty(new GeometryImpl<SRID>());
            return empty;
        }

        switch(ggRaw->getGeometryTypeId())
        {
        case gg::GEOS_POINT:
        {
            boost::shared_ptr<GeometryImpl<SRID> > pt(new PointImpl<SRID>());
            pt->pGeometry_ = boost::shared_ptr<gg::Point>(dynamic_cast<gg::Point*>(ggRaw), destroy);
            return pt;
        }
        case gg::GEOS_LINESTRING:
        {
            boost::shared_ptr<GeometryImpl<SRID> > ls(new LineStringImpl<SRID>());
            ls->pGeometry_ = boost::shared_ptr<gg::LineString>(dynamic_cast<gg::LineString*>(ggRaw), destroy);
            return ls;
        }
        case gg::GEOS_LINEARRING:
        {
            boost::shared_ptr<GeometryImpl<SRID> > ring(new LinearRingImpl<SRID>());
            ring->pGeometry_ = boost::shared_ptr<gg::LinearRing>(dynamic_cast<gg::LinearRing*>(ggRaw), destroy);
            return ring;
        }
        case gg::GEOS_POLYGON:
        {
            boost::shared_ptr<GeometryImpl<SRID> > poly(new PolygonImpl<SRID>());
            poly->pGeometry_ = boost::shared_ptr<gg::Polygon>(dynamic_cast<gg::Polygon*>(ggRaw), destroy);
            return poly;
        }
        case gg::GEOS_MULTIPOINT:
        {
            boost::shared_ptr<GeometryCollectionImpl<SRID> > mpt(new MultiPointImpl<SRID>());
            mpt->pGeometry_ = boost::shared_ptr<gg::MultiPoint>(dynamic_cast<gg::MultiPoint*>(ggRaw), destroy);
            return mpt;
        }
        case gg::GEOS_MULTILINESTRING:
        {
            boost::shared_ptr<GeometryCollectionImpl<SRID> > mls(new MultiLineStringImpl<SRID>());
            mls->pGeometry_ = boost::shared_ptr<gg::MultiLineString>(dynamic_cast<gg::MultiLineString*>(ggRaw), destroy);
            return mls;
        }
        case gg::GEOS_MULTIPOLYGON:
        {
            boost::shared_ptr<GeometryCollectionImpl<SRID> > mpoly(new MultiPolygonImpl<SRID>());
            mpoly->pGeometry_ = boost::shared_ptr<gg::MultiPolygon>(dynamic_cast<gg::MultiPolygon*>(ggRaw), destroy);
            return mpoly;
        }
        case gg::GEOS_GEOMETRYCOLLECTION:
        default:
        {
            boost::shared_ptr<GeometryCollectionImpl<SRID> > collection(new GeometryCollectionImpl<SRID>());
            collection->pGeometry_ = boost::shared_ptr<gg::GeometryCollection>(dynamic_cast<gg::GeometryCollection*>(ggRaw), destroy);
            return collection;
        }
        }
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::fromWKT(const std::string& wkt)
    {
        try {
            gg::Geometry* ggRaw = pWtkReader_->read(wkt);
            return fromRaw(ggRaw);
        } catch (geos::io::ParseException& pe) {
            return boost::shared_ptr<GeometryImpl<SRID> >(new GeometryImpl<SRID>());
        }
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createPoint(double x, double y)
    {
        boost::shared_ptr<GeometryImpl<SRID> > pt(new PointImpl<SRID>());
        pt->pGeometry_= boost::shared_ptr<gg::Point>(pFactory_->createPoint(gg::Coordinate(x, y))); //), destroyFunctor);
        return pt;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createPolygon(CoordinatePtrSequence const& coordinates)
    {
        boost::shared_ptr<GeometryImpl<SRID> > poly(new PolygonImpl<SRID>());
        gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
        BOOST_FOREACH(boost::shared_ptr<Coordinate> const c, coordinates)
        {
            seq->add(gg::Coordinate(c->getX(), c->getY()));
        }
        poly->pGeometry_ = boost::shared_ptr<gg::Polygon>(pFactory_->createPolygon(pFactory_->createLinearRing(seq), NULL));
        return poly;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createLineString(CoordinatePtrSequence const& coordinates)
    {
        boost::shared_ptr<GeometryImpl<SRID> > ls(new LineStringImpl<SRID>());
        gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
        BOOST_FOREACH(boost::shared_ptr<Coordinate> const c, coordinates)
        {
            seq->add(gg::Coordinate(c->getX(), c->getY()));
        }
        ls->pGeometry_ = boost::shared_ptr<gg::LineString>(pFactory_->createLineString(seq));
        return ls;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createLinearRing(CoordinatePtrSequence const& coordinates)
    {
        boost::shared_ptr<GeometryImpl<SRID> > ring(new LineStringImpl<SRID>());
        gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
        BOOST_FOREACH(boost::shared_ptr<Coordinate> const c, coordinates)
        {
            seq->add(gg::Coordinate(c->getX(), c->getY()));
        }
        ring->pGeometry_ = boost::shared_ptr<gg::LineString>(pFactory_->createLinearRing(seq));
        return ring;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createMultiPoint(CoordinatePtrSequence const& coordinates)
    {
        boost::shared_ptr<GeometryImpl<SRID> > mpt(new MultiPointImpl<SRID>());
        gg::CoordinateArraySequence seq;
        BOOST_FOREACH(boost::shared_ptr<Coordinate> const c, coordinates)
        {
            seq.add(gg::Coordinate(c->getX(), c->getY()));
        }
        mpt->pGeometry_ = boost::shared_ptr<gg::MultiPoint>(pFactory_->createMultiPoint(seq));
        return mpt;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createMultiLineString(std::vector<CoordinatePtrSequence> const& coordinates)
    {
        boost::shared_ptr<GeometryImpl<SRID> > mls(new MultiLineStringImpl<SRID>());
        std::vector<gg::Geometry*>* lsVec = new std::vector<gg::Geometry*>(coordinates.size());
        BOOST_FOREACH(CoordinatePtrSequence const& cseq, coordinates)
        {
            gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
            BOOST_FOREACH(CoordinatePtr const& coord, cseq)
            {
                seq->add(gg::Coordinate(coord->getX(), coord->getY()));
            }
            lsVec->push_back(pFactory_->createLineString(seq));
        }
        mls->pGeometry_ = boost::shared_ptr<gg::Geometry>(pFactory_->createMultiLineString(lsVec));
        return mls;
    }

    template <int SRID>
    boost::shared_ptr<GeometryImpl<SRID> > GeometryImpl<SRID>::createMultiPolygon(std::vector<CoordinatePtrSequence> const& coordinates)
    {
        boost::shared_ptr<GeometryImpl<SRID> > mpoly(new MultiPolygonImpl<SRID>());
        std::vector<gg::Geometry*>* polyVec = new std::vector<gg::Geometry*>(coordinates.size());
        BOOST_FOREACH(CoordinatePtrSequence const& cseq, coordinates)
        {
            gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
            BOOST_FOREACH(CoordinatePtr const& coord, cseq)
            {
                seq->add(gg::Coordinate(coord->getX(), coord->getY()));
            }
            polyVec->push_back(pFactory_->createPolygon(pFactory_->createLinearRing(seq), NULL));
        }
        mpoly->pGeometry_ = boost::shared_ptr<gg::Geometry>(pFactory_->createMultiPolygon(polyVec));
        return mpoly;
    }

} }
#endif // S3_GEOMETRY_IMPL_H_
