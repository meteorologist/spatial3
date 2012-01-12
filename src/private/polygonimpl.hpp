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

#ifndef S3_POLYGON_IMPL_H_
#define S3_POLYGON_IMPL_H_

// implementation
//
#include "geometryimpl.hpp"

// GEOC C++
//
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
// polygon related
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPolygon.h>

namespace gi = geos::io;
namespace gg = geos::geom;

namespace metno { namespace s3 {

    extern void destroy(gg::Geometry* g);

    // ## -------------------------------------------------------------------------------------------------------
    //     PolygonImpl
    // ## -------------------------------------------------------------------------------------------------------
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
    class MultiPolygonImpl : public GeometryCollectionImpl<SRID> {
    public:
        MultiPolygonImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createMultiPolygon());
        }

        boost::shared_ptr<GeometryImpl<SRID> > makeValid() const;
        boost::shared_ptr<GeometryImpl<SRID> > cleanGeometry() const;
    };

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

} }

#endif // S3_POLYGON_IMPL_H_
