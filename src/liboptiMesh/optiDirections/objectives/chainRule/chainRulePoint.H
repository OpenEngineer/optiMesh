#ifndef chainRulePoint_H
#define chainRulePoint_H

#include "point.H"
#include "tensorField.H"

#include "chainRuleFwd.H"

// Gradient coefficient meaning: if the tensor is right-multiplied by a vector, the resulting vector is equivalent to the gradient of the dot-product of the initial vector and that vector
//
//
// this & other:
//
//  this.x*other.x + this.y*other.y + this.z*other.z
//
// grad(this&other);
//
// [
//  this & dx(other) + other & dx(this);
//  this & dy(other) + other & dy(this);
//  this & dz(other) + other & dz(this)
// ]
//
// Can be split into two parts:
//
// [
//  dx(other);
//  dy(other);  * this + ...
//  dz(other);
// ]
//
// So components become:
//
// [
//   dcomp_x/dx, dcomp_y/dx, dcomp_z/dx;
//   dcomp_x/dy, dcomp_y/dy, dcomp_z/dy; 
//   dcomp_x/dz, dcomp_y/dz, dcomp_z/dz; 
// ]
//

namespace Foam {
namespace optiMesh {

class chainRulePoint 
{
  point value_;

  tensorField gradient_;

public:
  static label nPoints_;

  chainRulePoint();

  chainRulePoint(const point& value); // uses static member nPoints_

  chainRulePoint(const scalar& x, const scalar& y, const scalar& z); // uses static member nPoints_

  chainRulePoint(const point& value, const label& nPoints);

  chainRulePoint(const label& nPoints); // value=(0,0,0)

  chainRulePoint(const point& value, const label& id, const label& nPoints);

  chainRulePoint(const scalar& x, const scalar& y, const scalar& z, const label& id, 
      const label& nPoints);

  chainRulePoint(const point& value, const tensorField& gradient);

  ~chainRulePoint();

  vector value() const;

  tensor gradient(const label& id) const;

  label nPoints() const;


  crPoint& operator+=(const crPoint&);

  crPoint& operator+=(const point&);

  crPoint& operator-=(const crPoint&);

  crPoint& operator-=(const point&);

  crPoint& operator*=(const crScalar&);

  crPoint& operator*=(const scalar&);

  crPoint& operator/=(const crScalar&);

  crPoint& operator/=(const scalar&);


  friend crScalar operator&(const crPoint&, const crPoint&);

  friend crPoint operator^(const crPoint&, const crPoint&);
};

}
}

#endif
