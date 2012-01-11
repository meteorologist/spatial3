
#ifndef S3_GEOMETRY_H_
#define S3_GEOMETRY_H_

#include "coordinate.hpp"

// GEOC C++
//
#include <geos/geom/Geometry.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTWriter.h>

// boost
//
#include <boost/shared_ptr.hpp>

namespace gg = geos::geom;
namespace gi = geos::io;

namespace metno { namespace s3 {

    class LineStringD2;
    class GeometryFactory;

    typedef std::vector<Coordinate> CoordinateSequence;
    typedef std::vector<CoordinateSequence> CoordinateSequenceVector;

    class Geometry {

        friend class GeometryFactory;
        friend boost::shared_ptr<Geometry> fromRaw(gg::Geometry* raw);

    public:

        Geometry() {}

        Geometry(Geometry const &rhs)
        {
            pGeometry_ = rhs.pGeometry_;
            pFactory_  = rhs.pFactory_;
        }

        Geometry(Geometry const* g)
        {
            pGeometry_ = g->pGeometry_;
            pFactory_  = g->pFactory_;
        }

        Geometry &operator=(Geometry const &rhs)
        {
             if (this == &rhs)
                 return *this;

            pGeometry_ = rhs.pGeometry_;
            pFactory_  = rhs.pFactory_;
        }

        virtual ~Geometry();

        int getSRID() const;

        virtual bool isEmpty()     const;
        virtual bool isValid()     const;
        virtual bool isSimple()    const;
        virtual bool isRectangle() const;

        virtual double getArea() const;
        virtual double getLength() const;
        virtual Coordinate getCentroid() const;
        virtual Coordinate getLowerLeftCorner() const;
        virtual Coordinate getUpperRightCorner() const;
        virtual size_t getNumVertices() const;
        virtual size_t getNumGeometries() const;
        std::vector<Coordinate> getCoordinates() const;
        boost::shared_ptr<Geometry> getGeometryN(size_t n) const;

        // disjoint
        bool disjoint(Geometry const* g) const;
        bool disjoint(Geometry const& g) const;
        bool disjoint(boost::shared_ptr<Geometry> const& g) const;

        // touches
        bool touches(Geometry const* g) const;
        bool touches(Geometry const& g) const;
        bool touches(boost::shared_ptr<Geometry> const& g) const;

        // intersects
        bool intersects(Geometry const* g) const;
        bool intersects(Geometry const& g) const;
        bool intersects(boost::shared_ptr<Geometry> const& g) const;

        // crosses
        bool crosses(Geometry const* g) const;
        bool crosses(Geometry const& g) const;
        bool crosses(boost::shared_ptr<Geometry> const& g) const;

        // within
        bool within(Geometry const* g) const;
        bool within(Geometry const& g) const;
        bool within(boost::shared_ptr<Geometry> const& g) const;

        // contains
        bool contains(Geometry const* g) const;
        bool contains(Geometry const& g) const;
        bool contains(boost::shared_ptr<Geometry> const& g) const;

        // overlaps
        bool overlaps(Geometry const* g) const;
        bool overlaps(Geometry const& g) const;
        bool overlaps(boost::shared_ptr<Geometry> const& g) const;

        //equals
        bool equals(Geometry const* g) const;
        bool equals(Geometry const& g) const;
        bool equals(boost::shared_ptr<Geometry> const& g) const;

        // equalsExact
        bool equalsExact(Geometry const* g, double tolerance = 0) const;
        bool equalsExact(Geometry const& g, double tolerance = 0) const;
        bool equalsExact(boost::shared_ptr<Geometry> const& g, double tolerance = 0) const;

        // covers
        bool covers(Geometry const* g) const;
        bool covers(Geometry const& g) const;
        bool covers(boost::shared_ptr<Geometry> const& g) const;

        // coveredBy
        bool coveredBy(Geometry const* g) const;
        bool coveredBy(Geometry const& g) const;
        bool coveredBy(boost::shared_ptr<Geometry> const& g) const;


        // Intersection
        boost::shared_ptr<Geometry> Intersection(Geometry const* g) const;
        boost::shared_ptr<Geometry> Intersection(Geometry const& g) const;
        boost::shared_ptr<Geometry> Intersection(boost::shared_ptr<Geometry> const& g) const;

        // Union
        boost::shared_ptr<Geometry> Union(Geometry const* geometry) const;
        boost::shared_ptr<Geometry> Union(Geometry const& geometry) const;
        boost::shared_ptr<Geometry> Union(boost::shared_ptr<Geometry> const& pGeometry) const;

        // Difference
        virtual Geometry* Difference(Geometry const* geometry) const;
//        virtual boost::shared_ptr<Geometry> Difference(Geometry const& geometry) const;
//        virtual boost::shared_ptr<Geometry> Difference(boost::shared_ptr<Geometry> const& pGeometry) const;

        // SymDifference
//        boost::shared_ptr<Geometry> SymDifference(Geometry const* geometry) const;
//        boost::shared_ptr<Geometry> SymDifference(Geometry const& geometry) const;
//        boost::shared_ptr<Geometry> SymDifference(boost::shared_ptr<Geometry> const& pGeometry) const;

        virtual boost::shared_ptr<Geometry> Buffer(double distance) const;
        boost::shared_ptr<Geometry> Simplify(double tolerance) const;

        virtual std::string toWKT() const;

    protected:

//        virtual boost::shared_ptr<gg::Geometry> getGeosGeometry() { return pGeometry_; }

        GeometryFactory* pFactory_;
        boost::shared_ptr<gg::Geometry> pGeometry_;
    };

    class GeometryCollection : public Geometry {
    public:
        GeometryCollection();
    };

    class GeometryFactory {

        friend class Geometry;

    public:
        GeometryFactory();
        virtual ~GeometryFactory();

        virtual bool init() = 0;

        virtual boost::shared_ptr<Geometry> createPoint(double xAxis = 0.0, double yAxis = 0.0) = 0;
        virtual boost::shared_ptr<Geometry> createLineString(std::vector<Coordinate>& coordinates) = 0;
        virtual boost::shared_ptr<Geometry> createLinearRing(std::vector<Coordinate>& coordinates) = 0;
        virtual boost::shared_ptr<Geometry> createPolygon(std::vector<Coordinate>& coordinates) = 0;
//        virtual boost::shared_ptr<Geometry> createPolygon(boost::shared_ptr<LineString> const &ring) = 0;
        virtual boost::shared_ptr<Geometry> createMultiPoint(std::vector<Coordinate>& coordinates) = 0;
        virtual boost::shared_ptr<Geometry> createMultiLineString(CoordinateSequenceVector& coordinates) = 0;
        virtual boost::shared_ptr<Geometry> createMultiPolygon(CoordinateSequenceVector& coordinates) = 0;

        virtual boost::shared_ptr<Geometry> fromWKT(const std::string& wkt) = 0;
        virtual boost::shared_ptr<Geometry> fromRaw(boost::shared_ptr<Geometry>& raw) = 0;

    protected:

        virtual void destroyGeometry(Geometry* g);
        virtual void destroyGeometry(Geometry& g);
        virtual void destroyGeometry(boost::shared_ptr<Geometry>& g);

        virtual boost::shared_ptr<Geometry> createDefault() = 0;


        boost::shared_ptr<gg::PrecisionModel>  pm_;
        boost::shared_ptr<gg::GeometryFactory> geoFactory_;
        boost::shared_ptr<gi::WKTReader>       wktReader_;
        boost::shared_ptr<gi::WKTWriter>       wktWriter_;
    };

} }
#endif // S3GEOMETRY_H_
