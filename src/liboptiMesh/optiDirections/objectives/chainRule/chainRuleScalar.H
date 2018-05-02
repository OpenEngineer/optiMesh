#ifndef chainRuleScalar_H
#define chainRuleScalar_H

#include "label.H"
#include "scalar.H"
#include "point.H"
#include "vectorField.H"

#include "chainRuleFwd.H"

namespace Foam {
namespace optiMesh {

class chainRuleScalar
{
  friend chainRulePoint;

  scalar value_;

  vectorField gradient_;

public:
  static label nPoints_;

  chainRuleScalar(); // no initialization of member!

  chainRuleScalar(const scalar& value); // use the static nPoints_ member

  chainRuleScalar(const scalar& value, const label& nPoints);

  chainRuleScalar(const label& nPoints); // value == 0

  chainRuleScalar(const scalar& value, const vectorField& vf);

  ~chainRuleScalar();

  scalar value() const;

  vector gradient(const label& id) const;

  label nPoints() const;


  crScalar& operator+=(const crScalar&);

  crScalar& operator+=(const scalar&);

  crScalar& operator-=(const crScalar&);

  crScalar& operator-=(const scalar&);

  crScalar& operator*=(const crScalar&);

  crScalar& operator*=(const scalar&);

  crScalar& operator/=(const crScalar&);

  crScalar& operator/=(const scalar&);


  friend crScalar sqrt(const crScalar&);

  friend crScalar cbrt(const crScalar&);

  friend crScalar pow(const crScalar&, const scalar&);

  friend crScalar operator/(const scalar&, const crScalar&);
};

}
}

#endif
