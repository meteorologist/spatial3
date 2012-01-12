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

#ifndef S3_POINT_IMPL_H_
#define S3_POINT_IMPL_H_

// implementation
//
#include "geometryimpl.hpp"

// GEOC C++
//
#include <geos/geom/Geometry.h>
#include <geos/geom/GeometryFactory.h>
// point related
#include <geos/geom/Point.h>
#include <geos/geom/MultiPoint.h>

namespace gi = geos::io;
namespace gg = geos::geom;

namespace metno { namespace s3 {

    extern void destroy(gg::Geometry* g);

    // ## -------------------------------------------------------------------------------------------------------
    //     PointImpl
    // ## -------------------------------------------------------------------------------------------------------

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
    };

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
    //     MultiPointImpl
    // ## -------------------------------------------------------------------------------------------------------

    template <int SRID>
    class MultiPointImpl : public GeometryCollectionImpl<SRID> {
    public:
        MultiPointImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createMultiPoint());
        }
    };
}}
#endif // S3_POINT_IMPL_H_
