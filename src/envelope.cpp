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

#include "envelope.hpp"
#include "private/envelopeimpl.hpp"

namespace metno { namespace s3 {
Envelope::Envelope() : pImpl_(boost::shared_ptr<EnvelopeImpl>(new EnvelopeImpl()))
{
}

Envelope::Envelope(double x1, double x2, double y1, double y2)
    : pImpl_(boost::shared_ptr<EnvelopeImpl>(new EnvelopeImpl(x1, x2, y1, y2)))
{
}

Envelope::Envelope(const Coordinate &p1, const Coordinate &p2)
   : pImpl_(boost::shared_ptr<EnvelopeImpl>(new EnvelopeImpl(p1, p2)))
{
}

double Envelope::getArea() const
{
    return pImpl_->getArea();
}

double Envelope::getMaxX() const
{
    return pImpl_->getMaxX();
}

double Envelope::getMaxY() const
{
    return pImpl_->getMaxY();
}

double Envelope::getMinX() const
{
    return pImpl_->getMinX();
}

double Envelope::getMinY() const
{
    return pImpl_->getMinY();
}

bool Envelope::equals(const Envelope *other) const
{
    return pImpl_->equals(other->pImpl_.get());
}

std::string Envelope::toString() const
{
    return pImpl_->toString();
}

} }
