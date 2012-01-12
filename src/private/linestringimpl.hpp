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

#ifndef S3_LINESTRING_IMPL_H_
#define S3_LINESTRING_IMPL_H_

// implementation
//
#include "geometryimpl.hpp"

namespace gi = geos::io;
namespace gg = geos::geom;

namespace metno { namespace s3 {

    extern void destroy(gg::Geometry* g);

    // ## -------------------------------------------------------------------------------------------------------
    //     LineStringImpl
    // ## -------------------------------------------------------------------------------------------------------
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

        //        boost::shared_ptr<GeometryImpl> reverse();

    };

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
    //     MultiLineStringImpl
    // ## -------------------------------------------------------------------------------------------------------
    template <int SRID>
    class MultiLineStringImpl : public GeometryCollectionImpl<SRID> {
    public:
        MultiLineStringImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createMultiLineString());
        }
    };
} }

#endif // S3_LINESTRING_IMPL_H_
