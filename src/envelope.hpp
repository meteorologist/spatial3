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

#ifndef S3_ENVELOPE_H_
#define S3_ENVELOPE_H_

// implementation
//
#include "coordinate.hpp"

namespace metno { namespace s3 {

    class EnvelopeImpl;

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
        bool   equals(const Envelope *other) const;
        std::string toString() const;

    private:
        boost::shared_ptr<EnvelopeImpl> pImpl_;
    };

} }

#endif // S3_ENVELOPE_H_
