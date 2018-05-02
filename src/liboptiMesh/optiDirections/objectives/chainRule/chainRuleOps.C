#include "chainRuleScalar.H"
#include "chainRulePoint.H"
#include "chainRuleOps.H"


namespace Foam
{
namespace optiMesh 
{


tensor operator^(const tensor& a, const vector& b)
{
  return tensor(a.x() ^ b, a.y() ^ b, a.z() ^ b);
}


tensorField operator^(const tensorField& a, const vector& b)
{
  tensorField c(a.size());

  forAll(a, i) {
    c[i] = a[i] ^ b;
  }

  return c;
}


vectorField operator&(const tensorField& a, const vector& b)
{
  vectorField c(a.size());

  forAll(a, i) {
    c[i] = a[i] & b;
  }

  return c;
}


tensorField operator*(const vectorField& a, const vector& b)
{
  tensorField c(a.size());

  forAll(a, i) {
    c[i] = a[i] * b;
  }

  return c;
}


crScalar sqrt(const crScalar& a)
{
  scalar value = Foam::sqrt(a.value_);

  vectorField gradient = 0.5 * a.gradient_ / value;

  return crScalar(value, gradient);
}


scalar sqrt(const scalar& a)
{
  return Foam::sqrt(a);
}


crScalar cbrt(const crScalar& a)
{
  scalar value = Foam::cbrt(a.value_);

  vectorField gradient = (1.0/3.0) * a.gradient_ / (value * value);

  return crScalar(value, gradient);
}


scalar cbrt(const scalar& a)
{
  return Foam::cbrt(a);
}


crScalar mag(const crPoint& a)
{
  return sqrt(a&a);
}


scalar mag(const point& a)
{
  return Foam::mag(a);
}


crScalar pow(const crScalar& a, const scalar& p)
{
  if (a.value_ < 0.0) {
    FatalErrorInFunction << "input scalar is negative, check algorithm" << exit(FatalError);
  }

  scalar value1 = Foam::pow(a.value_, p-1); // more robust then pow(value_, p)

  vectorField gradient = p * a.gradient_ * value1;

  return crScalar(value1 * a.value_, gradient);
}


scalar pow(const scalar& a, const scalar& p) 
{
  if (a < 0.0) {
    FatalErrorInFunction << "input scalar is negative, check algorithm" << exit(FatalError);
  }

  return Foam::pow(a, p);
}


crScalar operator+(const crScalar& a, const crScalar& b)
{
  crScalar c(a);

  c += b;

  return c;
}


crScalar operator+(const crScalar& a, const scalar& b)
{
  crScalar c(a);

  c += b;

  return c;
}


crScalar operator+(const scalar& a, const crScalar& b)
{
  return b + a;
}


crScalar operator-(const crScalar& a, const crScalar& b)
{
  crScalar c(a);

  c -= b;

  return c;
}


crScalar operator-(const crScalar& a, const scalar& b)
{
  crScalar c(a);

  c -= b;

  return c;
}


crScalar operator-(const scalar& a, const crScalar& b)
{
  return ((b - a) *= -1.0);
}


crScalar operator*(const crScalar& a, const crScalar& b)
{
  crScalar c(a);

  c *= b;

  return c;
}


crScalar operator*(const crScalar& a, const scalar& b)
{
  crScalar c(a);

  c *= b;

  return c;
}


crScalar operator*(const scalar& a, const crScalar& b)
{
  return b*a;
}


crScalar operator/(const crScalar& a, const crScalar& b)
{
  crScalar c(a);

  c /= b;

  return c;
}


crScalar operator/(const crScalar& a, const scalar& b)
{
  crScalar c(a);

  c /= b;

  return c;
}


crScalar operator/(const scalar& a, const crScalar& b)
{
  scalar value = a / b.value_;

  vectorField gradient = -(value / b.value_) * b.gradient_;

  return crScalar(value, gradient);
}


crPoint operator+(const crPoint& a, const crPoint& b)
{
  crPoint c(a);

  c += b;

  return c;
}


crPoint operator-(const crPoint& a, const crPoint& b)
{
  crPoint c(a);

  c -= b;

  return c;
}


crPoint operator*(const crPoint& a, const crScalar& b)
{
  crPoint c(a);

  c *= b;

  return c;
}


crPoint operator*(const crScalar& a, const crPoint& b)
{
  return b*a;
}


crPoint operator*(const crPoint& a, const scalar& b)
{
  crPoint c(a);

  c *= b;

  return c;
}


crPoint operator*(const scalar& a, const crPoint& b)
{
  return b*a;
}


crPoint operator/(const crPoint& a, const crScalar& b)
{
  crPoint c(a);

  c /= b;

  return c;
}


crPoint operator/(const crPoint& a, const scalar& b)
{
  crPoint c(a);

  c /= b;

  return c;
}


crScalar operator&(const crPoint& a, const crPoint& b)
{
  scalar value = a.value_ & b.value_;

  vectorField gradient = (a.gradient_ & b.value_) +
                         (b.gradient_ & a.value_);

  return crScalar(value, gradient);
}


crPoint operator^(const crPoint& a, const crPoint& b)
{
  point value = a.value_ ^ b.value_;

  tensorField gradient = (a.gradient_ ^ b.value_) -
                         (b.gradient_ ^ a.value_);

  return crPoint(value, gradient);
}

}
}
