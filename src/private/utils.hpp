/*
 * spatial3
 *
 * (C) Copyright 2011, met.no
 *
 * Project Info:
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */

#ifndef S3_UTILS_IMPL_H_
#define S3_UTILS_IMPL_H_

// GEOS
//
#include <geos/geom/Geometry.h>

namespace gg = geos::geom;

namespace metno { namespace s3 {
    void destroy(gg::Geometry* g)
    {
        if(!g) return;

        const gg::GeometryFactory* factory = g->getFactory();
        factory->destroyGeometry(g);
    }

    // caller must take ownership of returned geometry
    gg::Geometry* add_point_linestring(gg::Geometry const* g, gg::Coordinate c)
    {
        gg::CoordinateSequence* coords = g->getCoordinates();
        coords->add(c);
        return g->getFactory()->createLineString(coords);
    }

    gg::Geometry* close2d_ring(gg::Geometry const* g)
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


    gg::Geometry* make_geos_friendly_ring(gg::Geometry const* r)
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

    gg::Geometry* make_geos_friendly_polygon(gg::Geometry* p)
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
    gg::Geometry* node_lines(gg::Geometry const* lines)
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

    gg::Geometry* extract_unique_points(gg::Geometry const* g)
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


    gg::Geometry* polygonize(const gg::Geometry* const* vgeoms, size_t ngeoms)
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

    gg::Geometry* build_area(gg::Geometry const* geom_in)
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

    gg::Geometry* make_valid_polygon(gg::Geometry const* p)
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

    gg::Geometry* make_clean_polygon(gg::Geometry const* poly_in)
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

    gg::Geometry* line_split_by_line(const gg::Geometry* line_in, const gg::Geometry* blade_in)
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

    gg::Geometry* polygon_split_by_line(const gg::Geometry* polygon, const gg::Geometry* line)
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

    gg::Geometry* collection_split_by_line(const gg::Geometry* coll_in, const gg::Geometry* blade_in)
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

} }

#endif // S3_UTILS_IMPL_H_
