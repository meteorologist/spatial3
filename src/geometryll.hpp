#ifndef S3_GEOMETRY_LL_H_
#define S3_GEOMETRY_LL_H_

// spatial3
//
#include "coordinate.hpp"

// boost
//
#include <boost/shared_ptr.hpp>

// STD
//
#include <vector>

namespace metno { namespace s3 {

    template <int SRID>
    class GeometryImpl;

    class EnvelopeImpl;

    const int SRID_WGS84 = 4362;

    class Envelope {
    public:
        Envelope();
        Envelope(double x1, double x2, double y1, double y2);
        Envelope(const Coordinate &p1, const Coordinate &p2);

        double getArea() const;
        double getMaxX() const;
        double getMaxY() const;
        double getMinX() const;
        double getMinY() const;
        bool equals(const Envelope *other) const;
        std::string toString() const;

    private:
        boost::shared_ptr<EnvelopeImpl> pImpl_;
    };

    class GeometryLL {
    public:
        GeometryLL();
        virtual ~GeometryLL();

        virtual std::string toWKT() const;

        virtual bool isEmpty()     const;
        virtual bool isValid()     const;
        virtual bool isRectangle() const;
        virtual bool isSimple()    const;

        virtual int          getSRID() const;
        virtual GeometryType getGeometryType() const;
        virtual size_t getNumVertices() const;
        virtual size_t getNumGeometries() const;
        virtual void getCoordinates(CoordinatePtrSequence& coordinates) const;
        virtual Coordinate getCentroid() const;
        Envelope getEnvelope() const;
        virtual Coordinate getLowerLeftCorner() const;
        virtual Coordinate getUpperRightCorner() const;
        virtual boost::shared_ptr<GeometryLL> getGeometryN(size_t n) const;
        DimensionType getDimension() const;

        virtual bool intersects(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool touches(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool crosses(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool contains(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool within(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool overlaps(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool covers(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool coveredBy(boost::shared_ptr<GeometryLL> const& g) const;
        virtual bool equals(boost::shared_ptr<GeometryLL> const& g) const;

        virtual boost::shared_ptr<GeometryLL> Intersection(boost::shared_ptr<GeometryLL> const& g) const;
        virtual boost::shared_ptr<GeometryLL> Union(boost::shared_ptr<GeometryLL> const g) const;
        virtual boost::shared_ptr<GeometryLL> Difference(boost::shared_ptr<GeometryLL> const g) const;
        virtual boost::shared_ptr<GeometryLL> Buffer(double distance) const;

        static boost::shared_ptr<GeometryLL> createPoint(double x = 0.0, double y = 0.0);
        static boost::shared_ptr<GeometryLL> createPolygon(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryLL> createLineString(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryLL> createLinearRing(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryLL> createMultiPoint(CoordinatePtrSequence const& coordinates);
        static boost::shared_ptr<GeometryLL> createMultiLineString(std::vector<CoordinatePtrSequence> const& coordinates);
        static boost::shared_ptr<GeometryLL> createMultiPolygon(std::vector<CoordinatePtrSequence> const& coordinates);
        static boost::shared_ptr<GeometryLL> fromWKT(const std::string& wkt);

        virtual void simplifyDP(double tolerance);
        virtual void simplifyTP(double tolerance);

        virtual boost::shared_ptr<GeometryLL> makeValid() const;
        virtual boost::shared_ptr<GeometryLL> cleanGeometry() const;
        boost::shared_ptr<GeometryLL> splitBy(boost::shared_ptr<GeometryLL> const& g) const;

    protected:
        boost::shared_ptr<GeometryLL> makeGeosFriendly() const;

        boost::shared_ptr<GeometryImpl<SRID_WGS84> > pImpl_;
    };


    class GeometryCollectionLL : public GeometryLL {

    public:
        GeometryCollectionLL();

    };

    class PointLL : public GeometryLL {

    public:

        PointLL();
        ~PointLL();

        double getLon();
        double getLat();
        void setLon(double lon);
        void setLat(double lat);
    };

    class LineStringLL : public GeometryLL {

    public:

        bool isRing()   const;
        bool isClosed() const;

        boost::shared_ptr<PointLL> getEndPoint() const;
        boost::shared_ptr<PointLL> getStartPoint() const;

//        boost::shared_ptr<Geometry> reverse();

    protected:

    };

    class LinearRingLL : public LineStringLL {
    public:
    };

    class PolygonLL : public GeometryLL {
    public:
        boost::shared_ptr<LineStringLL> getExteriorRing();
        boost::shared_ptr<GeometryLL> makeValid() const;
        boost::shared_ptr<GeometryLL> cleanGeometry() const;
    protected:

    };

    class MultiPointLL : public GeometryCollectionLL {
    public:
//        boost::shared_ptr<Coordinate> getCoordinateN(size_t n) const;
    protected:

    };

    class MultiLineStringLL : public GeometryCollectionLL {
    public:

    protected:

    };

    class MultiPolygonLL : public GeometryCollectionLL {
    public:
        boost::shared_ptr<GeometryLL> makeValid() const;
        boost::shared_ptr<GeometryLL> cleanGeometry() const;
    protected:

    };
} }

#endif // S3_GEOMETRY_LL_H_
