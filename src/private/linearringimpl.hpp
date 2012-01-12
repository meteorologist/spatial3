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

#ifndef S3_LINEARRING_IMPL_H_
#define S3_LINEARRING_IMPL_H_

// implementation
//
#include "geometryimpl.hpp"

namespace metno { namespace s3 {

    extern void destroy(gg::Geometry* g);

    // ## -------------------------------------------------------------------------------------------------------
    //     LinearRingImpl
    // ## -------------------------------------------------------------------------------------------------------
    template<int SRID>
    class LinearRingImpl : public LineStringImpl<SRID> {
    public:
        LinearRingImpl() {
            this->pGeometry_ = boost::shared_ptr<gg::Geometry>(this->pFactory_->createLinearRing());
        }
    };
} }

#endif // S3_LINEARRING_IMPL_H_
