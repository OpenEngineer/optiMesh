#include "chainRuleScalar.H"
#include "chainRulePoint.H"
#include "chainRuleOps.H"

using namespace Foam;

using namespace Foam::optiMesh;


label chainRulePoint::nPoints_ = 0;

chainRulePoint::chainRulePoint()
{}


chainRulePoint::chainRulePoint(const point& value, const label& nPoints) :
  chainRulePoint(value, -1, nPoints)
{
}


chainRulePoint::chainRulePoint(const point& value) :
  chainRulePoint(value, nPoints_)
{
  if (nPoints_ == 0) {
    FatalErrorInFunction <<
      "static member nPoints_ not intialized" << exit(FatalError);
  }
}


chainRulePoint::chainRulePoint(const scalar& x, const scalar& y, const scalar& z) :
  chainRulePoint(point(x,y,z), nPoints_)
{
  if (nPoints_ == 0) {
    FatalErrorInFunction <<
      "static member nPoints_ not intialized" << exit(FatalError);
  }
}


chainRulePoint::chainRulePoint(const label& nPoints) :
  chainRulePoint(point(0,0,0), nPoints)
{
}


chainRulePoint::chainRulePoint(const point& value, const label& id, const label& nPoints) :
  value_(value),
  gradient_(nPoints)
{
  forAll(gradient_, i) {
    gradient_[i] = (id == i) ? 
      tensor(1,0,0,  0,1,0,  0,0,1) :
      tensor(0,0,0,  0,0,0,  0,0,0);
  }
}


chainRulePoint::chainRulePoint(const scalar& x, const scalar& y, const scalar& z, 
    const label& id, const label& nPoints) :
  chainRulePoint(point(x,y,z), id, nPoints)
{}


chainRulePoint::chainRulePoint(const point& value, const tensorField& gradient) :
  value_(value),
  gradient_(gradient)
{}


chainRulePoint::~chainRulePoint()
{}


vector chainRulePoint::value() const
{
  return value_;
}


tensor chainRulePoint::gradient(const label& id) const
{
  return gradient_[id];
}


label chainRulePoint::nPoints() const
{
  return gradient_.size();
}


crPoint& chainRulePoint::operator+=(const crPoint& p)
{
  value_ += p.value_;

  gradient_ += p.gradient_;

  return *this;
}


crPoint& chainRulePoint::operator+=(const point& p)
{
  value_ += p;

  return *this;
}


crPoint& chainRulePoint::operator-=(const crPoint& p)
{
  value_ -= p.value_;

  gradient_ -= p.gradient_;

  return *this;
}


crPoint& chainRulePoint::operator-=(const point& p)
{
  value_ -= p;

  return *this;
}


crPoint& chainRulePoint::operator*=(const crScalar& s)
{
  gradient_ *= s.value_;

  gradient_ += (s.gradient_ * value_);

  value_ *= s.value_;

  return *this;
}


crPoint& chainRulePoint::operator*=(const scalar& s)
{
  value_ *= s;

  gradient_ *= s;

  return *this;
}


crPoint& chainRulePoint::operator/=(const crScalar& s)
{
  return (*this *= (1.0 / s));
}


crPoint& chainRulePoint::operator/=(const scalar& s)
{
  return (*this *= (1.0 / s));
}
