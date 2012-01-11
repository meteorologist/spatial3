#include "coordinate.hpp"

// GEOS
//
#include "geos/operation/valid/IsValidOp.h"
#include "geos/simplify/DouglasPeuckerSimplifier.h"

namespace metno { namespace s3 {

void Coordinate::setX(double x)
{
    coordinate_.x = x;
}

double Coordinate::getX() const
{
    return coordinate_.x;
}

void Coordinate::setY(double y)
{
    coordinate_.y = y;
}

double Coordinate::getY() const
{
    return coordinate_.y;
}

std::string Coordinate::toString() const
{
    return coordinate_.toString();
}

bool Coordinate::equals(const Coordinate& coord) const
{
    return coordinate_.equals(coord.coordinate_);
}


double Coordinate::distance(const Coordinate& coord) const
{
    return coordinate_.distance(coord.coordinate_);
}

bool Coordinate::isValid() const
{
    geos::geom::Coordinate c(getX(), getY());
    return geos::operation::valid::IsValidOp::isValid(c);
}


bool CoordinatePtrSequence::hasRepeatedPoints() const
{
    gg::CoordinateArraySequence seq;
    for(size_t i = 0; i < size(); ++i)
        seq.add(gg::Coordinate(at(i)->getX(), at(i)->getY()));
    return seq.hasRepeatedPoints();
}

void CoordinatePtrSequence::removeRepeatedPoints()
{
    gg::CoordinateArraySequence seq;
    for(size_t i = 0; i < size(); ++i)
        seq.add(gg::Coordinate(at(i)->getX(), at(i)->getY()));

    gg::CoordinateSequence& cleared(seq.removeRepeatedPoints());

    clear();

    for(size_t i = 0; i < cleared.getSize(); ++i)
        push_back(CoordinatePtr(new Coordinate(cleared.getX(i), cleared.getY(i))));
}

void CoordinatePtrSequence::reverse()
{
    gg::CoordinateArraySequence seq;
    for(size_t i = 0; i < size(); ++i)
        seq.add(gg::Coordinate(at(i)->getX(), at(i)->getY()));

    gg::CoordinateSequence::reverse(&seq);

    clear();

    for(size_t i = 0; i < seq.getSize(); ++i)
        push_back(CoordinatePtr(new Coordinate(seq.getX(i), seq.getY(i))));
}

void CoordinatePtrSequence::simplify(double distance)
{
    assert(0);
    gg::CoordinateArraySequence seq;
    for(size_t i = 0; i < size(); ++i)
        seq.add(gg::Coordinate(at(i)->getX(), at(i)->getY()));



    for(size_t i = 0; i < seq.getSize(); ++i)
        push_back(CoordinatePtr(new Coordinate(seq.getX(i), seq.getY(i))));
}

std::string CoordinatePtrSequence::toString() const
{
    gg::CoordinateArraySequence seq;
    for(size_t i = 0; i < size(); ++i)
        seq.add(gg::Coordinate(at(i)->getX(), at(i)->getY()));
    return seq.toString();
}

// friend free operators
//
bool operator== (Coordinate const& lhs, Coordinate const& rhs)
{
    return lhs.equals(rhs);
}

bool operator!= (Coordinate const& lhs, Coordinate const& rhs)
{
    return !(lhs == rhs);
}

bool operator<  (Coordinate const& lhs, Coordinate const& rhs)
{
    gg::CoordinateLessThen clt;
    return clt(lhs.coordinate_, rhs.coordinate_);
}

bool operator>  (Coordinate const& lhs, Coordinate const& rhs)
{
    return !(lhs < rhs);
}

}} // namespaces
