#include "geometry.hpp"

// spatialite
//
#include "backend/headers/spatialite/sqlite3.h"
#include "backend/headers/spatialite/gaiageo.h"
#include "backend/headers/spatialite.h"

// GEOS
//
#include <geos/geom/Envelope.h>
#include <geos/simplify/DouglasPeuckerSimplifier.h>

// BOOST
//
#include <boost/foreach.hpp>
#include <boost/pointer_cast.hpp>

// STD
//
#include <limits>

namespace metno { namespace s3 {

    static void release(gg::Geometry* g)
    {
        if(!g)
            return;

        //        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "
        //                 <<g->toString()
        //                 <<std::endl;

        const gg::GeometryFactory* factory = g->getFactory();
        factory->destroyGeometry(g);
        //       std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" PASS "<<std::endl;
    }

// ## -------------------------------------------------------------------------------------------------------
//     Geometry
// ## -------------------------------------------------------------------------------------------------------

    Geometry::~Geometry()
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" PASS "<<std::endl;
    }

    int Geometry::getSRID() const
    {
        if(pGeometry_.get())
            return pGeometry_->getSRID();
        else
            return -1;
    }

    bool Geometry::isEmpty() const
    {
        if(pGeometry_.get())
            return pGeometry_->isEmpty();
        else
            return false;
    }

    bool Geometry::isValid() const
    {
        if(pGeometry_.get())
            return pGeometry_->isValid();
        else
            return -1;
    }

    bool Geometry::isSimple() const
    {
        if(pGeometry_.get())
            return pGeometry_->isSimple();
        else
            return false;
    }

    bool Geometry::isRectangle() const
    {
        if(pGeometry_.get())
            return pGeometry_->isRectangle();
        else
            return false;
    }

    double Geometry::getArea() const
    {
        if(pGeometry_.get())
            return pGeometry_->getArea();
        else
            return -1;
    }

    double Geometry::getLength() const
    {
        if(pGeometry_.get())
            return pGeometry_->getLength();
        else
            return -1;
    }

    Coordinate Geometry::getCentroid() const
    {
        if(pGeometry_.get()) {
            gg::Coordinate ggC;
            pGeometry_->getCentroid(ggC);
            return Coordinate(ggC.x, ggC.y);
        } else {
            return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
    }

    Coordinate Geometry::getUpperRightCorner() const
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

    Coordinate Geometry::getLowerLeftCorner() const
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
        if(pGeometry_.get()) {
//            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            if(pGeometry_->getEnvelopeInternal() != 0) {
//                std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" ENVELOPE: "<<pGeometry_->getEnvelopeInternal()->toString()<<std::endl;
                return Coordinate(pGeometry_->getEnvelopeInternal()->getMinX(), pGeometry_->getEnvelopeInternal()->getMinY());
            }
        }

        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
        return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }

    size_t Geometry::getNumVertices() const
    {
        if(pGeometry_.get())
            return pGeometry_->getNumPoints();
        else
            return 0;
    }

    size_t Geometry::getNumGeometries() const
    {
        if(pGeometry_.get())
            return pGeometry_->getNumGeometries();
        else
            return 0;
    }

    std::vector<Coordinate> Geometry::getCoordinates() const
    {
        std::vector<Coordinate> coordinates;
        if(pGeometry_.get() == 0 || pGeometry_->isEmpty())
            return coordinates;

        gg::CoordinateSequence* seq = pGeometry_->getCoordinates();
        std::vector<gg::Coordinate> coords;
        seq->toVector(coords);
        BOOST_FOREACH(gg::Coordinate const& c, coords)
        {
            coordinates.push_back(Coordinate(c.x, c.y));
        }
        delete seq;
        return coordinates;
    }

    std::string Geometry::toWKT() const
    {
        if(pGeometry_.get())
            return pGeometry_->toString();
        else
            return std::string();
    }

    bool Geometry::disjoint(const Geometry *g) const
    {
        if(pGeometry_.get())
            return pGeometry_->disjoint(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::disjoint(Geometry const& g) const
    {
        return disjoint(&g);
    }

    bool Geometry::disjoint(boost::shared_ptr<Geometry> const& g) const
    {
        return disjoint(g.get());
    }

    bool Geometry::touches(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->touches(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::touches(Geometry const& g) const
    {
        return touches(&g);
    }

    bool Geometry::touches(boost::shared_ptr<Geometry> const& g) const
    {
        return touches(g.get());
    }

    bool Geometry::intersects(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->intersects(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::intersects(Geometry const& g) const
    {
        return intersects(&g);
    }

    bool Geometry::intersects(boost::shared_ptr<Geometry> const& g) const
    {
        return intersects(g.get());
    }

    bool Geometry::crosses(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->crosses(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::crosses(Geometry const& g) const
    {
        return crosses(&g);
    }

    bool Geometry::crosses(boost::shared_ptr<Geometry> const& g) const
    {
        return crosses(g.get());
    }

    bool Geometry::within(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->within(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::within(Geometry const& g) const
    {
        return within(&g);
    }

    bool Geometry::within(boost::shared_ptr<Geometry> const& g) const
    {
        return within(g.get());
    }

    bool Geometry::contains(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->contains(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::contains(Geometry const& g) const
    {
        return contains(&g);
    }

    bool Geometry::contains(boost::shared_ptr<Geometry> const& g) const
    {
        return contains(g.get());
    }

    bool Geometry::overlaps(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->overlaps(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::overlaps(Geometry const& g) const
    {
        return overlaps(&g);
    }

    bool Geometry::overlaps(boost::shared_ptr<Geometry> const& g) const
    {
        return overlaps(g.get());
    }

    bool Geometry::equals(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->equals(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::equals(Geometry const& g) const
    {
        return equals(&g);
    }

    bool Geometry::equals(boost::shared_ptr<Geometry> const& g) const
    {
        return equals(g.get());
    }

    bool Geometry::covers(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->covers(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::covers(Geometry const& g) const
    {
        return covers(&g);
    }

    bool Geometry::covers(boost::shared_ptr<Geometry> const& g) const
    {
        return covers(g.get());
    }

    bool Geometry::coveredBy(Geometry const* g) const
    {
        if(pGeometry_.get())
            return pGeometry_->coveredBy(g->pGeometry_.get());
        else
            return false;
    }

    bool Geometry::coveredBy(Geometry const& g) const
    {
        return coveredBy(&g);
    }

    bool Geometry::coveredBy(boost::shared_ptr<Geometry> const& g) const
    {
        return coveredBy(g.get());
    }

    bool Geometry::equalsExact(Geometry const* g, double tolerance) const
    {
        if(pGeometry_.get())
            return pGeometry_->equalsExact(g->pGeometry_.get(), tolerance);
        else
            return false;
    }

    bool Geometry::equalsExact(Geometry const& g, double tolerance) const
    {
        return equalsExact(&g, tolerance);
    }

    bool Geometry::equalsExact(boost::shared_ptr<Geometry> const& g, double tolerance) const
    {
        return equalsExact(g.get(), tolerance);
    }

    boost::shared_ptr<Geometry> Geometry::getGeometryN(size_t n) const
    {
        if(pGeometry_.get())
            assert(0);//pGeometry_->getGeometryN(n);
        else
            return boost::shared_ptr<Geometry>();
    }

    boost::shared_ptr<Geometry> Geometry::Intersection(Geometry const* geometry) const
    {
        std::auto_ptr<gg::Geometry> ggGeom(pGeometry_->intersection(geometry->pGeometry_.get()));
        boost::shared_ptr<Geometry> raw = pFactory_->createDefault();
        raw->pGeometry_ = boost::shared_ptr<gg::Geometry>(ggGeom->clone());
        return pFactory_->fromRaw(raw);
    }

    boost::shared_ptr<Geometry> Geometry::Intersection(Geometry const& geometry) const
    {
        return Intersection(&geometry);
    }

    boost::shared_ptr<Geometry> Geometry::Intersection(boost::shared_ptr<Geometry> const& geometry) const
    {
        return Intersection(geometry.get());
    }

    boost::shared_ptr<Geometry> Geometry::Union(Geometry const* geometry) const
    {
        std::auto_ptr<gg::Geometry> ggGeom(pGeometry_->Union(geometry->pGeometry_.get()));
        boost::shared_ptr<Geometry> raw = pFactory_->createDefault();
        raw->pGeometry_ = boost::shared_ptr<gg::Geometry>(ggGeom->clone());
        return pFactory_->fromRaw(raw);
    }

    boost::shared_ptr<Geometry> Geometry::Union(Geometry const& geometry) const
    {
        return Union(&geometry);
    }

    boost::shared_ptr<Geometry> Geometry::Union(boost::shared_ptr<Geometry> const& geometry) const
    {
        return Union(geometry.get());
    }

    Geometry* Geometry::Difference(Geometry const* geometry) const
    {
        assert(geometry);
        assert(0);
        return 0;
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;

//        gg::Geometry* geos = pGeometry_->difference(geometry->pGeometry_.get());

//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;

//        return pFactory_->fromRaw(geosgeometry);
    }

//    boost::shared_ptr<Geometry> Geometry::Difference(Geometry const& geometry) const
//    {
////        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;
//        return Difference(&geometry);
//    }

//    boost::shared_ptr<Geometry> Geometry::Difference(boost::shared_ptr<Geometry> const& geometry) const
//    {
////        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;
//        return Difference(geometry.get());
//    }

//    boost::shared_ptr<Geometry> Geometry::SymDifference(Geometry const* geometry) const
//    {
//        std::auto_ptr<gg::Geometry> ggGeom(pGeometry_->symDifference(geometry->pGeometry_.get()));
//        boost::shared_ptr<Geometry> raw = pFactory_->createDefault();
//        raw->pGeometry_ = boost::shared_ptr<gg::Geometry>(ggGeom->clone());
//        return pFactory_->fromRaw(raw);
//    }

//    boost::shared_ptr<Geometry> Geometry::SymDifference(Geometry const& geometry) const
//    {
//        return SymDifference(&geometry);
//    }

//    boost::shared_ptr<Geometry> Geometry::SymDifference(boost::shared_ptr<Geometry> const& geometry) const
//    {
//        return SymDifference(geometry.get());
//    }

    boost::shared_ptr<Geometry> Geometry::Buffer(double distance) const
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;
        std::auto_ptr<gg::Geometry> ggGeom(pGeometry_->buffer(distance));
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;
        assert(pFactory_);
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;
        boost::shared_ptr<Geometry> raw(new Geometry());
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : CHECK "<<std::endl;
        raw->pGeometry_ = boost::shared_ptr<gg::Geometry>(ggGeom->clone());
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : WKT : "<<raw->pGeometry_->toString()<<std::endl;
        boost::shared_ptr<Geometry> tmp = pFactory_->fromRaw(raw);
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<std::endl;
        return tmp;
    }

    boost::shared_ptr<Geometry> Geometry::Simplify(double tolerance) const
    {
        std::auto_ptr<gg::Geometry> ggGeom
                = geos::simplify::DouglasPeuckerSimplifier::simplify(pGeometry_.get(), tolerance);

        boost::shared_ptr<Geometry> raw = pFactory_->createDefault();
        raw->pGeometry_ = boost::shared_ptr<gg::Geometry>(ggGeom->clone());
        return pFactory_->fromRaw(raw);
    }

// ## -------------------------------------------------------------------------------------------------------
//     GeometryCollection
// ## -------------------------------------------------------------------------------------------------------

    GeometryCollection::GeometryCollection() { }

// ## -------------------------------------------------------------------------------------------------------
//     GeometryFactory
// ## -------------------------------------------------------------------------------------------------------

    GeometryFactory::GeometryFactory() { }

    GeometryFactory::~GeometryFactory() { }

    void GeometryFactory::destroyGeometry(Geometry* g)
    {
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" PASS "<<std::endl;
        if(g && geoFactory_.get())
            geoFactory_->destroyGeometry(g->pGeometry_.get());
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" PASS "<<std::endl;
    }

    void GeometryFactory::destroyGeometry(Geometry& g)
    {
        destroyGeometry(&g);
    }

    void GeometryFactory::destroyGeometry(boost::shared_ptr<Geometry>& g)
    {
        destroyGeometry(g.get());
    }
} }
