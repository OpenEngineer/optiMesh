#ifndef chainRuleOps_H
#define chainRuleOps_H

#include "tensor.H"
#include "point.H"
#include "scalar.H"

#include "chainRuleFwd.H"


namespace Foam {
namespace optiMesh {

// special operator
tensor operator^(const tensor&, const vector&);
tensorField operator^(const tensorField&, const vector&);
vectorField operator&(const tensorField&, const vector&);
tensorField operator*(const vectorField&, const vector&);

// unary operators
crScalar sqrt(const crScalar&);
scalar sqrt(const scalar&); // just to have it in same namespace
crScalar cbrt(const crScalar&);
scalar cbrt(const scalar&); // just to have it in same namespace
crScalar mag(const crPoint&);
scalar mag(const point&); // just to have it in same namespace
crScalar pow(const crScalar&, const scalar&);
scalar pow(const scalar&, const scalar&);
// just use a&a for magSq


// binary operators
crScalar operator+(const crScalar&, const crScalar&);
crScalar operator+(const crScalar&, const scalar&);
crScalar operator+(const scalar&, const crScalar&);
crScalar operator-(const crScalar&, const crScalar&);
crScalar operator-(const crScalar&, const scalar&);
crScalar operator-(const scalar&,   const crScalar&);
crScalar operator*(const crScalar&, const crScalar&);
crScalar operator*(const crScalar&, const scalar&);
crScalar operator*(const scalar&,   const crScalar&);
crScalar operator/(const crScalar&, const crScalar&);
crScalar operator/(const scalar&,   const crScalar&);
crScalar operator/(const crScalar&, const scalar&);

crPoint  operator+(const crPoint&,  const crPoint&);
crPoint  operator-(const crPoint&,  const crPoint&);
crPoint  operator*(const crPoint&,  const crScalar&);
crPoint  operator*(const crScalar&, const crPoint&);
crPoint  operator*(const crPoint&,  const scalar&);
crPoint  operator*(const scalar&,   const crPoint&);
crPoint  operator/(const crPoint&,  const crScalar&);
crPoint  operator/(const crPoint&,  const scalar&);

crScalar operator&(const crPoint&,  const crPoint&);
crPoint  operator^(const crPoint&,  const crPoint&);

}
}

#endif
