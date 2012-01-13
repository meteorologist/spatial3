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

#ifndef S3_GEOMETRY_IMPL_H_
#define S3_GEOMETRY_IMPL_H_

#include "coordinate.hpp"
#include "envelopeimpl.hpp"

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

namespace gi = geos::io;
namespace gg = geos::geom;
namespace go = geos::operation;

namespace metno { namespace s3 {

    template <int SRID>
    class PointImpl;

    template <int SRID>
    class MultiPointImpl;

    template <int SRID>
    class LineStringImpl;

    template <int SRID>
    class MultiLineStringImpl;

    template <int SRID>
    class LinearRingImpl;

    template <int SRID>
    class PolygonImpl;

    template <int SRID>
    class MultiPolygonImpl;

    extern void destroy(gg::Geometry* g);

    // caller must take ownership of returned geometry
    extern gg::Geometry* make_geos_friendly_polygon(gg::Geometry* p);
    extern gg::Geometry* polygon_split_by_line(const gg::Geometry* polygon, const gg::Geometry* line);
    extern gg::Geometry* collection_split_by_line(const gg::Geometry* coll_in, const gg::Geometry* blade_in);

// ## -------------------------------------------------------------------------------------------------------
//     GeometryImpl
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
