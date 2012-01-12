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

#ifndef S3_ENVELOPE_IMPL_H_
#define S3_ENVELOPE_IMPL_H_



// implementation
//
#include "coordinate.hpp"

// GEOS
//
#include <geos/geom/Envelope.h>

namespace gg = geos::geom;

namespace metno { namespace s3 {

    // ## -------------------------------------------------------------------------------------------------------
    //     EnvelopeImpl
    // ## -------------------------------------------------------------------------------------------------------

    class EnvelopeImpl {
    public:
        EnvelopeImpl()
            : envelope_(new gg::Envelope())
        {}
        EnvelopeImpl(double x1, double x2, double y1, double y2)
            : envelope_(new gg::Envelope(x1, x2, y1, y2))
        {}

        EnvelopeImpl(const Coordinate &p1, const Coordinate &p2)
            : envelope_(new gg::Envelope(gg::Coordinate(p1.getX(), p1.getY()), gg::Coordinate(p2.getX(), p2.getY())))
        {}

        double getArea() const {
            return envelope_->getArea();
        }

        double getMaxX() const {
            return envelope_->getMaxX();
        }

        double getMaxY() const {
            return envelope_->getMaxY();
        }
        double getMinX() const {
            return envelope_->getMinX();
        }

        double getMinY() const {
            return envelope_->getMinY();
        }

        std::string toString() const {
            return envelope_->toString();
        }

        bool centre(Coordinate &centre) const {
            gg::Coordinate ggC;
            bool result = envelope_->centre(ggC);
            centre.setX(ggC.x);
            centre.setY(ggC.y);
            return result;
        }

        void setToNull() {
            envelope_->setToNull();
        }

        bool isNull() const {
            return envelope_->isNull();
        }

        bool equals(const EnvelopeImpl *other) const {
            return envelope_->equals(other->envelope_.get());
        }

        double getWidth() const {
            return envelope_->getWidth();
        }

        double getHeight() const {
            return envelope_->getHeight();
        }

        bool intersection(const EnvelopeImpl &env, EnvelopeImpl &result) const {
            return envelope_->intersection(*env.envelope_, *result.envelope_);
        }

        bool intersects(EnvelopeImpl const* other) const {
            return envelope_->intersects(other->envelope_.get());
        }

        bool intersects(Coordinate const& c) const {
            return envelope_->intersects(gg::Coordinate(c.getX(), c.getY()));
        }

    protected:
    private:
        boost::shared_ptr<gg::Envelope> envelope_;
    };
} }

#endif // S3_ENVELOPE_IMPL_H_
