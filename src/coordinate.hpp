#ifndef S3_COORDINATE_H_
#define S3_COORDINATE_H_

// GEOS
//
#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/index/quadtree/Quadtree.h>

// BOOST
//
#include <boost/shared_ptr.hpp>

// STD
//
#include <iostream>

namespace gg = geos::geom;

namespace metno { namespace s3 {

    class Coordinate;
    typedef boost::shared_ptr<Coordinate> CoordinatePtr;

    class Coordinate {
    public:

        Coordinate(double x = 0.0, double y = 0.0) : coordinate_(x, y){ }
        Coordinate(Coordinate const& c) : coordinate_(c.getX(), c.getY()){ }
        Coordinate(Coordinate const* c) : coordinate_(c->getX(), c->getY()){ }
        Coordinate(boost::shared_ptr<Coordinate> const& c) : coordinate_(c->getX(), c->getY()){ }

        virtual ~Coordinate() {}

        virtual bool isValid() const;

        virtual void setX(double x);
        virtual double getX() const;
        virtual void setY(double y);
        virtual double getY() const;
        virtual double distance(const Coordinate& coord) const;
        virtual bool equals(const Coordinate& coord) const;

        virtual std::string toString() const;

        friend bool operator== (Coordinate const& lhs, Coordinate const& rhs);
        friend bool operator!= (Coordinate const& lhs, Coordinate const& rhs);
        friend bool operator<  (Coordinate const& lhs, Coordinate const& rhs);
        friend bool operator>  (Coordinate const& lhs, Coordinate const& rhs);


    protected:
        gg::Coordinate coordinate_;
    };

    class CoordinatePtrSequence : public std::vector<CoordinatePtr> {
    public:
        bool hasRepeatedPoints() const;
        void removeRepeatedPoints();
        void reverse();
        void simplify(double distance);

        std::string toString() const;
    };

    enum GeometryType {
        // none
        S3_EMPTY,
        // a point
        S3_POINT,
        // a linestring
        S3_LINESTRING,
        // a linear ring (linestring with 1st point == last point)
        S3_LINEARRING,
        // a polygon
        S3_POLYGON,
        // a collection of points
        S3_MULTIPOINT,
        // a collection of linestrings
        S3_MULTILINESTRING,
        // a collection of polygons
        S3_MULTIPOLYGON,
        // a collection of heterogeneus geometries
        S3_GEOMETRYCOLLECTION
    };

    enum DimensionType {
        S3_DIM_POINT,
        S3_DIM_LINE,
        S3_DIM_SURFACE
    };

    class GeometryLL;
    class Envelope;
    class TedLocation;

    class Quad3Index {
    public:
        Quad3Index();
        ~Quad3Index();

        void insert(boost::shared_ptr<GeometryLL> geometry);
        void query(const boost::shared_ptr<Envelope> searchEnv, std::vector<void*>& ret);

    private:
        boost::shared_ptr<geos::index::quadtree::Quadtree> quad3_;
    };

} }

#endif
