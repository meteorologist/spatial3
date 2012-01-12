#ifndef S3_REGION_TEMPLATE_H_
#define S3_REGION_TEMPLATE_H_

// BOOST
//
#include <boost/shared_ptr.hpp>

namespace gg = geos::geom;

namespace metno { namespace s3 {

    template <typename GeometryType>
    class Region {
    public:
        explicit Region() { }
        explicit Region(std::string const& name, int id) : name_(name), id_(id) { }

        ~Region() {}

        void clear();

        bool        isValid() const;
        int         getId() const;
        void        setId(int id);
	void        setName(std::string const& name);
        std::string getName() const;
        double      getLength() const;
        double      getArea() const;
	size_t      getNumVertices() const;


        Coordinate getCentroid() const;
        Coordinate getLowerLeftCorner() const;
        Coordinate getUpperRightCorner() const;

        bool isInside(Coordinate const& c, int tolerance = 0) const;
        bool isInside(Region<GeometryType> const& r, int tolerance = 0) const;
        bool isIdentical(Region<GeometryType> const& r, int treshold = 0) const;

//        virtual void setCorners(const std::vector<Coordinate> &c);

    protected:
        int id_;
        std::string name_;
        boost::shared_ptr<GeometryType> pGeometry_;
    };

    template <typename GeometryType>
    void Region<GeometryType>::clear() {
        if(pGeometry_.get())
            pGeometry_  = boost::shared_ptr<GeometryType>(new GeometryType());
    }

    template <typename GeometryType>
    bool Region<GeometryType>::isValid() const
    {
        if(pGeometry_.get())
            return pGeometry_->isValid();
        else
            return false;
    }

    template <typename GeometryType>
    int Region<GeometryType>::getId() const { return id_; }

    template <typename GeometryType>
    void Region<GeometryType>::setId(int id)
    {
//        std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" : set id = "<<id_<<std::endl;
        id_ = id;
    }
    
    template <typename GeometryType>
    std::string Region<GeometryType>::getName() const { return name_; }

    template <typename GeometryType>
    void Region<GeometryType>::setName(std::string const& name) { name_ = name; }
    
    template <typename GeometryType>
    double Region<GeometryType>::getArea() const
    {
        if(pGeometry_.get())
            return pGeometry_->getArea();
        else
            return 0.0;
    }

    template <typename GeometryType>
    double Region<GeometryType>::getLength() const
    {
        if(pGeometry_.get())
            return pGeometry_->getLength();
        else
            return 0.0;
    }

    template <typename GeometryType>
    size_t Region<GeometryType>::getNumVertices() const
    {
        if(pGeometry_.get())
            return pGeometry_->getNumVertices();
        else
            return 0;
    }
    
    template <typename GeometryType>
    Coordinate Region<GeometryType>::getCentroid() const
    {
        if(pGeometry_.get())
            return pGeometry_->getCentroid();
        else
            return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    }

    template <typename GeometryType>
    Coordinate Region<GeometryType>::getUpperRightCorner() const
    {
        if(pGeometry_.get()) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            return pGeometry_->getUpperRightCorner();
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
    }

    template <typename GeometryType>
    Coordinate Region<GeometryType>::getLowerLeftCorner() const
    {
        if(pGeometry_.get()) {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            return pGeometry_->getLowerLeftCorner();
        } else {
            std::cerr<<"["<<__FILE__<<"]"<<__FUNCTION__<<"@"<<__LINE__<<" CHECK "<<std::endl;
            return Coordinate(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
    }

    template <typename GeometryType>
    bool Region<GeometryType>::isInside(Coordinate const& c, int tolerance) const
    {
        if(pGeometry_.get()) {
            boost::shared_ptr<GeometryType> buffer = boost::dynamic_pointer_cast<GeometryType>(pGeometry_->Buffer(tolerance));
            boost::shared_ptr<GeometryType> point  = boost::dynamic_pointer_cast<GeometryType>(pGeometry_->getFactory()->createPoint(c.getX(), c.getY()));
            return buffer->coveredBy(point);
        }
        return false;
    }

    template <typename GeometryType>
    bool Region<GeometryType>::isInside(Region<GeometryType> const& r, int tolerance) const
    {
        if(pGeometry_.get()) {
            boost::shared_ptr<GeometryType> buffer = boost::dynamic_pointer_cast<GeometryType>(pGeometry_->Buffer(tolerance));
            return buffer->coveredBy(r.pGeometry_);
        }
        return false;
    }

    template <typename GeometryType>
    bool Region<GeometryType>::isIdentical(Region<GeometryType> const& lhs, int threshold) const
    {
        if(!lhs.isValid() || !isValid())
            return false;
        if(pGeometry_.get()) {
            return pGeometry_->equals(lhs.pGeometry_);
        }

        return false;
    }
} }

#endif

