#include "chainRuleScalar.H"


using namespace Foam;

using namespace Foam::optiMesh;


label chainRuleScalar::nPoints_ = 0;


chainRuleScalar::chainRuleScalar()
{}


chainRuleScalar::chainRuleScalar(const scalar& value, const label& nPoints) :
  value_(value),
  gradient_(nPoints)
{
  forAll(gradient_, i) {
    gradient_[i] = vector(0,0,0);
  }
}


chainRuleScalar::chainRuleScalar(const label& nPoints) :
  chainRuleScalar(0.0, nPoints)
{
}


chainRuleScalar::chainRuleScalar(const scalar& value) :
  chainRuleScalar(value, nPoints_)
{
  if (nPoints_ == 0) {
    FatalErrorInFunction <<
      "static member nPoints_ not yet initialized" << exit(FatalError);
  }
}


chainRuleScalar::chainRuleScalar(const scalar& value, const vectorField& gradient) :
  value_(value),
  gradient_(gradient)
{}


chainRuleScalar::~chainRuleScalar()
{}


scalar chainRuleScalar::value() const
{
  return value_;
}


vector chainRuleScalar::gradient(const label& id) const
{
  return gradient_[id];
}


label chainRuleScalar::nPoints() const
{
  return gradient_.size();
}


crScalar& chainRuleScalar::operator+=(const crScalar& s)
{
  value_ += s.value_;

  gradient_ += s.gradient_;

  return *this;
}


crScalar& chainRuleScalar::operator+=(const scalar& s)
{
  value_ += s;

  return *this;
}


crScalar& chainRuleScalar::operator-=(const crScalar& s)
{
  value_ -= s.value_;

  gradient_ -= s.gradient_;

  return *this;
}


crScalar& chainRuleScalar::operator-=(const scalar& s)
{
  value_ -= s;

  return *this;
}


crScalar& chainRuleScalar::operator*=(const crScalar& s)
{
  gradient_ *= s.value_;

  gradient_ += value_ * s.gradient_;

  value_ *= s.value_;

  return *this;
}


crScalar& chainRuleScalar::operator*=(const scalar& s)
{
  value_ *= s;

  gradient_ *= s;

  return *this;
}


crScalar& chainRuleScalar::operator/=(const crScalar& s)
{
  scalar newValue = value_ / s.value_;

  gradient_ -= s.gradient_ * newValue;

  gradient_ /= s.value_;

  value_ = newValue;

  return *this;
}


crScalar& chainRuleScalar::operator/=(const scalar& s)
{
  return (*this *= (1.0/s));
}
