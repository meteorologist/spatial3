/**
 * \file PolarStereographic.hpp
 * \brief Header for GeographicLib::PolarStereographic class
 *
 * Copyright (c) Charles Karney (2008, 2009, 2010, 2011) <charles@karney.com>
 * and licensed under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 **********************************************************************/

#if !defined(GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP)
#define GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP \
  "$Id: 5cb57eb4c444c8a0f691f887fce1e17e1ed3dca6 $"

#include <GeographicLib/Constants.hpp>

namespace GeographicLib {

  /**
   * \brief Polar Stereographic Projection
   *
   * Implementation taken from the report,
   * - J. P. Snyder,
   *   <a href="http://pubs.er.usgs.gov/usgspubs/pp/pp1395"> Map Projections: A
   *   Working Manual</a>, USGS Professional Paper 1395 (1987),
   *   pp. 160&ndash;163.
   *
   * This is a straightforward implementation of the equations in Snyder except
   * that Newton's method is used to invert the projection.
   **********************************************************************/
  class GEOGRAPHIC_EXPORT PolarStereographic {
  private:
    typedef Math::real real;
    // _Cx used to be _C but g++ 3.4 has a macro of that name
    real _a, _f, _e2, _e, _e2m, _Cx, _c;
    real _k0;
    static const real tol_;
    static const real overflow_;
    static const int numit_ = 5;
    // Return e * atanh(e * x) for f >= 0, else return
    // - sqrt(-e2) * atan( sqrt(-e2) * x) for f < 0
    inline real eatanhe(real x) const throw() {
      return _f >= 0 ? _e * Math::atanh(_e * x) : - _e * std::atan(_e * x);
    }
  public:

    /**
     * Constructor for a ellipsoid with
     *
     * @param[in] a equatorial radius (meters).
     * @param[in] f flattening of ellipsoid.  Setting \e f = 0 gives a sphere.
     *   Negative \e f gives a prolate ellipsoid.  If \e f > 1, set flattening
     *   to 1/\e f.
     * @param[in] k0 central scale factor.
     *
     * An exception is thrown if either of the axes of the ellipsoid is
     * not positive \e a or if \e k0 is not positive.
     **********************************************************************/
    PolarStereographic(real a, real f, real k0);

    /**
     * Set the scale for the projection.
     *
     * @param[in] lat (degrees) assuming \e northp = true.
     * @param[in] k scale at latitude \e lat (default 1).
     *
     * This allows a "latitude of true scale" to be specified.  An exception is
     * thrown if \e k is not positive or if \e lat is not in the range (-90,
     * 90].
     **********************************************************************/
    void SetScale(real lat, real k = real(1));

    /**
     * Forward projection, from geographic to polar stereographic.
     *
     * @param[in] northp the pole which is the center of projection (true means
     *   north, false means south).
     * @param[in] lat latitude of point (degrees).
     * @param[in] lon longitude of point (degrees).
     * @param[out] x easting of point (meters).
     * @param[out] y northing of point (meters).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added.  \e lat should be in the range
     * (-90, 90] for \e northp = true and in the range [-90, 90) for \e northp
     * = false; \e lon should be in the range [-180, 360].
     **********************************************************************/
    void Forward(bool northp, real lat, real lon,
                 real& x, real& y, real& gamma, real& k) const throw();

    /**
     * Reverse projection, from polar stereographic to geographic.
     *
     * @param[in] northp the pole which is the center of projection (true means
     *   north, false means south).
     * @param[in] x easting of point (meters).
     * @param[in] y northing of point (meters).
     * @param[out] lat latitude of point (degrees).
     * @param[out] lon longitude of point (degrees).
     * @param[out] gamma meridian convergence at point (degrees).
     * @param[out] k scale of projection at point.
     *
     * No false easting or northing is added.  The value of \e lon returned is
     * in the range [-180, 180).
     **********************************************************************/
    void Reverse(bool northp, real x, real y,
                 real& lat, real& lon, real& gamma, real& k) const throw();

    /**
     * PolarStereographic::Forward without returning the convergence and scale.
     **********************************************************************/
    void Forward(bool northp, real lat, real lon,
                 real& x, real& y) const throw() {
      real gamma, k;
      Forward(northp, lat, lon, x, y, gamma, k);
    }

    /**
     * PolarStereographic::Reverse without returning the convergence and scale.
     **********************************************************************/
    void Reverse(bool northp, real x, real y,
                 real& lat, real& lon) const throw() {
      real gamma, k;
      Reverse(northp, x, y, lat, lon, gamma, k);
    }

    /** \name Inspector functions
     **********************************************************************/
    ///@{
    /**
     * @return \e a the equatorial radius of the ellipsoid (meters).  This is
     *   the value used in the constructor.
     **********************************************************************/
    Math::real MajorRadius() const throw() { return _a; }

    /**
     * @return \e f the flattening of the ellipsoid.  This is the value used in
     *   the constructor.
     **********************************************************************/
    Math::real Flattening() const throw() { return _f; }

    /**
     * <b>DEPRECATED</b>
     * @return \e r the inverse flattening of the ellipsoid.
     **********************************************************************/
    Math::real InverseFlattening() const throw() { return 1/_f; }

    /**
     * The central scale for the projection.  This is the value of \e k0 used
     * in the constructor and is the scale at the pole unless overridden by
     * PolarStereographic::SetScale.
     **********************************************************************/
    Math::real CentralScale() const throw() { return _k0; }
    ///@}

    /**
     * A global instantiation of PolarStereographic with the WGS84 ellipsoid
     * and the UPS scale factor.  However, unlike UPS, no false easting or
     * northing is added.
     **********************************************************************/
    static const PolarStereographic UPS;
  };

} // namespace GeographicLib

#endif  // GEOGRAPHICLIB_POLARSTEREOGRAPHIC_HPP
