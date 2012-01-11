#include "geometryll.hpp"
#include "geometryimpl.hpp"

// SPATIALITE AMALGAMATION
//
#include <spatialite/sqlite3.h>
#include <spatialite/gaiageo.h>
#include <spatialite.h>

// GraphicLib - temporary helpers
//
#include <GeographicLib/PolygonArea.hpp>

// GEOS
//
#include <geos/geom/Point.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>

// BOOST
//
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/lexical_cast.hpp>

// STD
//
#include <string>
#include <vector>
#include <limits>

namespace metno { namespace s3 {

// ## -------------------------------------------------------------------------------------------------------
//     GeometryLL
// ## -------------------------------------------------------------------------------------------------------

GeometryLL::GeometryLL()
    : pImpl_(boost::shared_ptr<GeometryImpl<SRID_WGS84> >())
{
}

GeometryLL::~GeometryLL()
{
}

std::string GeometryLL::toWKT() const
{
    return pImpl_->toWKT();
}

int GeometryLL::getSRID() const
{
    return pImpl_->getSRID();
}

GeometryType GeometryLL::getGeometryType() const
{
    return pImpl_->getGeometryType();
}

bool GeometryLL::isEmpty() const
{
    return pImpl_->isEmpty();
}

bool GeometryLL::isValid() const
{
    return pImpl_->isValid();
}

bool GeometryLL::isRectangle() const
{
    return pImpl_->isRectangle();
}

bool GeometryLL::isSimple() const
{
    return pImpl_->isSimple();
}

bool GeometryLL::equals(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->equals(g->pImpl_);
}

size_t GeometryLL::getNumVertices() const
{
    return pImpl_->getNumVertices();
}

size_t GeometryLL::getNumGeometries() const
{
    return pImpl_->getNumGeometries();
}

void GeometryLL::getCoordinates(CoordinatePtrSequence& coordinates) const
{
    return pImpl_->getCoordinates(coordinates);
}

Coordinate GeometryLL::getCentroid() const
{
    return pImpl_->getCentroid();
}

Envelope GeometryLL::getEnvelope() const
{
    return Envelope(pImpl_->getEnvelope().getMinX(),
                    pImpl_->getEnvelope().getMaxX(),
                    pImpl_->getEnvelope().getMinY(),
                    pImpl_->getEnvelope().getMaxY());
}

Coordinate GeometryLL::getLowerLeftCorner() const
{
    return pImpl_->getLowerLeftCorner();
}

Coordinate GeometryLL::getUpperRightCorner() const
{
    return pImpl_->getUpperRightCorner();
}

bool GeometryLL::contains(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->contains(g->pImpl_);
}

bool GeometryLL::intersects(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->intersects(g->pImpl_);
}

bool GeometryLL::touches(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->touches(g->pImpl_);
}

bool GeometryLL::crosses(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->crosses(g->pImpl_);
}

bool GeometryLL::within(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->within(g->pImpl_);
}

bool GeometryLL::overlaps(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->overlaps(g->pImpl_);
}

bool GeometryLL::covers(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->covers(g->pImpl_);
}

bool GeometryLL::coveredBy(boost::shared_ptr<GeometryLL> const& g) const
{
    return pImpl_->coveredBy(g->pImpl_);
}

boost::shared_ptr<GeometryLL> GeometryLL::getGeometryN(size_t n) const
{
    boost::shared_ptr<GeometryLL> result(new GeometryLL());
    result->pImpl_ = pImpl_->getGeometryN(n);
    return result;
}

DimensionType GeometryLL::getDimension() const
{
    return pImpl_->getDimension();
}

boost::shared_ptr<GeometryLL> GeometryLL::Intersection(boost::shared_ptr<GeometryLL> const& g) const
{
    boost::shared_ptr<GeometryLL> result(new GeometryLL());
    result->pImpl_ = pImpl_->Intersection(g->pImpl_);
    return result;
}

boost::shared_ptr<GeometryLL> GeometryLL::Union(boost::shared_ptr<GeometryLL> const g) const
{
    boost::shared_ptr<GeometryLL> result(new GeometryLL());
    result->pImpl_ = pImpl_->Union(g->pImpl_);
    return result;
}

boost::shared_ptr<GeometryLL> GeometryLL::Difference(boost::shared_ptr<GeometryLL> const g) const
{
    boost::shared_ptr<GeometryLL> result(new GeometryLL());
    result->pImpl_ = pImpl_->Difference(g->pImpl_);
    return result;
}

boost::shared_ptr<GeometryLL> GeometryLL::Buffer(double distance) const
{
    boost::shared_ptr<GeometryLL> result(new GeometryLL());
    result->pImpl_ = pImpl_->Buffer(distance);
    return result;
}

boost::shared_ptr<GeometryLL> GeometryLL::createPoint(double x, double y)
{
    boost::shared_ptr<GeometryLL > pt(new PointLL());
    pt->pImpl_ = GeometryImpl<SRID_WGS84>::createPoint(x, y); //), destroyFunctor);
    return pt;
}

boost::shared_ptr<GeometryLL> GeometryLL::createPolygon(CoordinatePtrSequence const& coordinates)
{
    boost::shared_ptr<GeometryLL> poly(new PolygonLL());
    poly->pImpl_ = GeometryImpl<SRID_WGS84>::createPolygon(coordinates);
    return poly;
}

boost::shared_ptr<GeometryLL> GeometryLL::createLineString(CoordinatePtrSequence const& coordinates)
{
    boost::shared_ptr<GeometryLL> ls(new LineStringLL());
    ls->pImpl_ = GeometryImpl<SRID_WGS84>::createLineString(coordinates);
    return ls;
}

boost::shared_ptr<GeometryLL> GeometryLL::createLinearRing(CoordinatePtrSequence const& coordinates)
{
    boost::shared_ptr<GeometryLL> ring(new LinearRingLL());
    ring->pImpl_ = GeometryImpl<SRID_WGS84>::createLinearRing(coordinates);
    return ring;
}

boost::shared_ptr<GeometryLL> GeometryLL::createMultiPoint(CoordinatePtrSequence const& coordinates)
{
    boost::shared_ptr<GeometryLL> mpt(new MultiPointLL());
    mpt->pImpl_ = GeometryImpl<SRID_WGS84>::createMultiPoint(coordinates);
    return mpt;
}

boost::shared_ptr<GeometryLL> GeometryLL::createMultiLineString(std::vector<CoordinatePtrSequence> const& coordinates)
{
    boost::shared_ptr<GeometryLL> mls(new MultiLineStringLL());
    mls->pImpl_ = GeometryImpl<SRID_WGS84>::createMultiLineString(coordinates);
    return mls;
}

boost::shared_ptr<GeometryLL> GeometryLL::createMultiPolygon(std::vector<CoordinatePtrSequence> const& coordinates)
{
    boost::shared_ptr<GeometryLL> mpoly(new MultiPolygonLL());
    mpoly->pImpl_ = GeometryImpl<SRID_WGS84>::createMultiPolygon(coordinates);
    return mpoly;
}

boost::shared_ptr<GeometryLL> GeometryLL::fromWKT(const std::string& wkt)
{
    boost::shared_ptr<GeometryLL> poly(new PolygonLL());
    poly->pImpl_ = GeometryImpl<SRID_WGS84>::fromWKT(wkt);
    return poly;
}

//   double GeometryLL::getLength() const
//   {
//       std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT - START "<<std::endl;
//       GeometryFactoryLL* factory = boost::dynamic_pointer_cast<GeometryFactoryLL>(pFactory_);
//       double length = factory->getGodesicLengthFor(*this);
//       std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT - START "<<std::endl;
//       return length;
//   }

//   double GeometryLL::getArea() const
//   {
//       std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT - START "<<std::endl;
//       if(pGeometry_.get() == 0 || pGeometry_->isEmpty() || !pGeometry_->isValid() || pGeometry_->getArea() == 0)
//           return 0.0;
//       GeographicLib::PolygonArea polyarea(GeographicLib::Geodesic::WGS84);
//       double area = 0.0;
//       double perimeter = 0.0;
//       polyarea.Compute(false, true, perimeter, area);
//       std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT - END "<<std::endl;
//       return area;
//   }

void GeometryLL::simplifyDP(double tolerance)
{
    pImpl_->simplifyDP(tolerance);
}

void GeometryLL::simplifyTP(double tolerance)
{
    pImpl_->simplifyTP(tolerance);
}

boost::shared_ptr<GeometryLL> GeometryLL::makeGeosFriendly() const
{
    boost::shared_ptr<GeometryLL> geom;

    switch(getGeometryType()) {
    case S3_POLYGON:
        geom = boost::shared_ptr<GeometryLL>(new PolygonLL());
        break;
    case S3_MULTIPOLYGON:
        geom = boost::shared_ptr<GeometryLL>(new MultiPolygonLL());
        break;
    default:
        geom = boost::shared_ptr<GeometryLL>(new GeometryLL());
        break;
    }

    geom->pImpl_ = this->pImpl_->makeGeosFriendly();

    return geom;
}

boost::shared_ptr<GeometryLL> GeometryLL::makeValid() const
{
    boost::shared_ptr<GeometryLL> geom(new GeometryLL());
    geom->pImpl_ = this->pImpl_->makeValid();
    return geom;
}

boost::shared_ptr<GeometryLL> GeometryLL::cleanGeometry() const
{
    boost::shared_ptr<GeometryLL> geom(new GeometryLL());
    geom->pImpl_ = this->pImpl_->cleanGeometry();
    return geom;
}

boost::shared_ptr<GeometryLL> GeometryLL::splitBy(boost::shared_ptr<GeometryLL> const& g) const
{
    boost::shared_ptr<GeometryLL> geom(new GeometryLL());
    geom->pImpl_ = this->pImpl_->splitBy(g->pImpl_);
    return geom;
}

// ## -------------------------------------------------------------------------------------------------------
//     Envelope
// ## -------------------------------------------------------------------------------------------------------
Envelope::Envelope() : pImpl_(boost::shared_ptr<EnvelopeImpl>(new EnvelopeImpl())) { }

Envelope::Envelope(double x1, double x2, double y1, double y2)
    : pImpl_(boost::shared_ptr<EnvelopeImpl>(new EnvelopeImpl(x1, x2, y1, y2))) { }

Envelope::Envelope(const Coordinate &p1, const Coordinate &p2)
    : pImpl_(boost::shared_ptr<EnvelopeImpl>(new EnvelopeImpl(p1, p2))) {}

double Envelope::getArea() const {
    return pImpl_->getArea();
}

double Envelope::getMaxX() const {
    return pImpl_->getMaxX();
}

double Envelope::getMaxY() const {
    return pImpl_->getMaxY();
}

double Envelope::getMinX() const {
    return pImpl_->getMinX();
}

double Envelope::getMinY() const {
    return pImpl_->getMinY();
}

bool Envelope::equals(const Envelope *other) const {
    return pImpl_->equals(other->pImpl_.get());
}

std::string Envelope::toString() const {
    return pImpl_->toString();
}

// ## -------------------------------------------------------------------------------------------------------
//     PointLL
// ## -------------------------------------------------------------------------------------------------------

PointLL::PointLL()
{
//    std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT "<<std::endl;
    pImpl_ = boost::shared_ptr<GeometryImpl<SRID_WGS84> >(new PointImpl<SRID_WGS84>());
//    std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT "<<pImpl_->toWKT()<<std::endl;
}

PointLL::~PointLL()
{
}

double PointLL::getLon()
{
    boost::shared_ptr<PointImpl<SRID_WGS84> > pt = boost::dynamic_pointer_cast<PointImpl<SRID_WGS84> >(pImpl_);
    return pt->getX();
}

double PointLL::getLat()
{
    PointImpl<SRID_WGS84>* pt = dynamic_cast<PointImpl<SRID_WGS84>* >(pImpl_.get());
    return pt->getY();
}

void PointLL::setLon(double newlon)
{
    PointImpl<SRID_WGS84>* pt = dynamic_cast<PointImpl<SRID_WGS84>* >(pImpl_.get());
    pt->setX(newlon);
}

void PointLL::setLat(double newlat)
{
    boost::shared_ptr<PointImpl<SRID_WGS84> > pt = boost::dynamic_pointer_cast<PointImpl<SRID_WGS84> >(pImpl_);
    pt->setY(newlat);
}

// ## -------------------------------------------------------------------------------------------------------
//     SpatialDBInMemory
// ## -------------------------------------------------------------------------------------------------------

    class SpatialDBInMemory {
    public:
        SpatialDBInMemory(std::string const& dbname);
        ~SpatialDBInMemory();

        bool query(std::string const& query, std::vector<std::vector<std::string> >& results);

    protected:
        std::string  dbname_;
        sqlite3      *handle_;
        sqlite3_stmt *stmt_;

    };

    bool SpatialDBInMemory::query(std::string const& query, std::vector<std::vector<std::string> >& results)
    {
        sqlite3_stmt *statement;

        if(sqlite3_prepare_v2(handle_, query.c_str(), -1, &statement, 0) == SQLITE_OK) {
            int cols = sqlite3_column_count(statement);
            int result = 0;
            while(true) {
                result = sqlite3_step(statement);

                if(result == SQLITE_ROW) {
                    std::vector<std::string> values;
                    for(int col = 0; col < cols; col++) {
                        char * ptr = (char*)sqlite3_column_text(statement, col);
                        std::string val = ptr ? std::string(ptr) : std::string();
                        values.push_back(val);
                    }
                    results.push_back(values);
                } else {
                    break;
                }
            }

            sqlite3_finalize(statement);

            return true;
        } else {
            std::string error = sqlite3_errmsg(handle_);
            if(error != "not an error") std::cerr<<query<<" "<<error<<std::endl;
            results.clear();
            return false;
        }
    }

    SpatialDBInMemory::SpatialDBInMemory(std::string const& dbname)
        : dbname_(dbname)
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT "<<std::endl;
        spatialite_init(0);
        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "
                <<" SQLite version = "<<sqlite3_libversion ()<<"    "
                <<" SpatiaLite version = "<<spatialite_version()
                <<std::endl;

        if(sqlite3_open_v2 (":memory:", &handle_, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL) != SQLITE_OK) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "
                     <<" cannot open in memory db: "<<" error : "<<sqlite3_errmsg(handle_)<<std::endl;
            sqlite3_close(handle_);
            assert(0);
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "
                     <<" opened in memory db "<<std::endl;
        }

        std::vector<std::vector<std::string> > results;
        if(query("SELECT InitSpatialMetadata()", results)) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "
                     <<" InitSpatialMetadata() -> OK "<<std::endl;
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : "
                     <<" InitSpatialMetadata() -> NOK "<<std::endl;
            sqlite3_close(handle_);
            assert(0);
        }

    }

    SpatialDBInMemory::~SpatialDBInMemory()
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT "<<std::endl;
        sqlite3_close(handle_);
        spatialite_cleanup();
        sqlite3_reset_auto_extension();
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK POINT "<<std::endl;
    }

//   boost::shared_ptr<Geometry> GeometryFactoryLL::createLinearRing(std::vector<Coordinate>& coordinates)
//   {
//       boost::shared_ptr<LineStringLL> ls(new LineStringLL());
//       gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
//       BOOST_FOREACH(Coordinate const coord, coordinates)
//       {
//           seq->add(gg::Coordinate(coord.getX(), coord.getY()));
//       }
//       ls->pFactory_ = this;
//       ls->pGeometry_ = boost::shared_ptr<gg::Geometry>(geoFactory_->createLinearRing(seq));
//       return ls;
//   }


//   boost::shared_ptr<Geometry> GeometryFactoryLL::createMultiPoint(std::vector<Coordinate>& coordinates)
//   {
//       boost::shared_ptr<MultiPointLL> poly(new MultiPointLL());
//       gg::CoordinateArraySequence seq(coordinates.size());
//       BOOST_FOREACH(Coordinate const cd2, coordinates)
//       {
//           seq.add(gg::Coordinate(cd2.getX(), cd2.getY()));
//       }
//       poly->pFactory_ = this;
//       poly->pGeometry_ = boost::shared_ptr<gg::MultiPoint>(geoFactory_->createMultiPoint(seq));
//       return poly;
//   }

//   boost::shared_ptr<Geometry> GeometryFactoryLL::createMultiLineString(CoordinateSequenceVector& coordinates)
//   {
//       boost::shared_ptr<MultiLineStringLL> mls(new MultiLineStringLL());
//       std::vector<gg::Geometry*>* lsVec = new std::vector<gg::Geometry*>(coordinates.size());
//       BOOST_FOREACH(CoordinateSequence const& cseq, coordinates)
//       {
//           gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
//           BOOST_FOREACH(Coordinate const& coord, cseq)
//           {
//               seq->add(gg::Coordinate(coord.getX(), coord.getY()));
//           }
//           lsVec->push_back(geoFactory_->createLineString(seq));
//       }
//       mls->pFactory_ = this;
//       mls->pGeometry_ = boost::shared_ptr<gg::Geometry>(geoFactory_->createMultiLineString(lsVec));
//       return mls;
//   }

//   boost::shared_ptr<Geometry> GeometryFactoryLL::createMultiPolygon(CoordinateSequenceVector& coordinates)
//   {
//       boost::shared_ptr<MultiPolygonLL> mpoly(new MultiPolygonLL());
//       std::vector<gg::Geometry*>* polyVec = new std::vector<gg::Geometry*>(coordinates.size());
//       BOOST_FOREACH(CoordinateSequence const& cseq, coordinates)
//       {
//           gg::CoordinateSequence *seq = new gg::CoordinateArraySequence();
//           BOOST_FOREACH(Coordinate const& coord, cseq)
//           {
//               seq->add(gg::Coordinate(coord.getX(), coord.getY()));
//           }
//           polyVec->push_back(geoFactory_->createPolygon(geoFactory_->createLinearRing(seq), NULL));
//       }
//       mpoly->pFactory_ = this;
//       mpoly->pGeometry_ = boost::shared_ptr<gg::Geometry>(geoFactory_->createMultiPolygon(polyVec));
//       return mpoly;
//   }

// ## -------------------------------------------------------------------------------------------------------
//     LineStringLL
// ## -------------------------------------------------------------------------------------------------------

   bool LineStringLL::isRing() const
   {
       LineStringImpl<SRID_WGS84>* ls = dynamic_cast<LineStringImpl<SRID_WGS84>* >(pImpl_.get());
       return ls->isRing();
   }

   bool LineStringLL::isClosed() const
   {
       LineStringImpl<SRID_WGS84>* ls = dynamic_cast<LineStringImpl<SRID_WGS84>* >(pImpl_.get());
       return ls->isClosed();
   }

   boost::shared_ptr<PointLL> LineStringLL::getEndPoint() const
   {
       LineStringImpl<SRID_WGS84>* ls = dynamic_cast<LineStringImpl<SRID_WGS84>* >(pImpl_.get());
       boost::shared_ptr<PointImpl<SRID_WGS84> > end = ls->getEndPoint();
       return boost::dynamic_pointer_cast<PointLL>(GeometryLL::createPoint(end->getX(), end->getY()));
   }

   boost::shared_ptr<PointLL> LineStringLL::getStartPoint() const
   {
       LineStringImpl<SRID_WGS84>* ls = dynamic_cast<LineStringImpl<SRID_WGS84>* >(pImpl_.get());
       boost::shared_ptr<PointImpl<SRID_WGS84> > start = ls->getStartPoint();
       return boost::dynamic_pointer_cast<PointLL>(GeometryLL::createPoint(start->getX(), start->getY()));
   }

//   boost::shared_ptr<Geometry> LineStringLL::reverse()
//   {
//       boost::shared_ptr<gg::LineString> ls = boost::dynamic_pointer_cast<gg::LineString>(pGeometry_);

//       gg::Geometry* geom = ls->reverse();
//       gg::CoordinateSequence* seq = geom->getCoordinates();
//       std::vector<gg::Coordinate> coordinates;
//       seq->toVector(coordinates);

//       std::vector<Coordinate> coords;
//       BOOST_FOREACH(gg::Coordinate const coord, coordinates)
//       {
//           coords.push_back(Coordinate(coord.x, coord.y));
//       }
//       pGeometry_->getFactory()->destroyGeometry(geom);
//       return pFactory_->createLineString(coords);
//   }

// ## -------------------------------------------------------------------------------------------------------
//     PolygonLL
// ## -------------------------------------------------------------------------------------------------------

    boost::shared_ptr<LineStringLL> PolygonLL::getExteriorRing()
    {
        PolygonImpl<SRID_WGS84>* poly = dynamic_cast<PolygonImpl<SRID_WGS84>* >(pImpl_.get());
        boost::shared_ptr<GeometryImpl<SRID_WGS84> > shell = poly->getExteriorRing();
        if(shell.get()) {
            CoordinatePtrSequence seq;
            shell->getCoordinates(seq);
            return boost::dynamic_pointer_cast<LineStringLL>(GeometryLL::createLineString(seq));
        } else {
            return boost::shared_ptr<LineStringLL>(new LineStringLL());
        }
    }

    boost::shared_ptr<GeometryLL> PolygonLL::makeValid() const
    {
        boost::shared_ptr<PolygonLL> geom(new PolygonLL());
        geom->pImpl_ = this->pImpl_->makeValid();
        return geom;
    }

    boost::shared_ptr<GeometryLL> PolygonLL::cleanGeometry() const
    {
        boost::shared_ptr<PolygonLL> geom(new PolygonLL());
        geom->pImpl_ = this->pImpl_->cleanGeometry();
        return geom;
    }

// ## -------------------------------------------------------------------------------------------------------
//     MultiPolygonLL
// ## -------------------------------------------------------------------------------------------------------

    boost::shared_ptr<GeometryLL> MultiPolygonLL::makeValid() const
    {
        boost::shared_ptr<MultiPolygonLL> geom(new MultiPolygonLL());
        geom->pImpl_ = this->pImpl_->makeValid();
        return geom;
    }

    boost::shared_ptr<GeometryLL> MultiPolygonLL::cleanGeometry() const
    {
        boost::shared_ptr<MultiPolygonLL> geom(new MultiPolygonLL());
        geom->pImpl_ = this->pImpl_->cleanGeometry();
        return geom;
    }

// ## -------------------------------------------------------------------------------------------------------
//     GeometryCollectionLL
// ## -------------------------------------------------------------------------------------------------------

    GeometryCollectionLL::GeometryCollectionLL() : GeometryLL() { }

} }
