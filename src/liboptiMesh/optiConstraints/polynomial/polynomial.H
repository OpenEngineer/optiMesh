#ifndef polynomial_H
#define polynomial_H

#include "optiConstraint.H"

// system/smoothMeshDict entry:
//
// {
//   type polynomial;
//   set p0;
//   origin (0 0 0);
//   axis (0 0 1);  // project point onto the 
//   n (0 1 0); // t = n ^ axis, keep t coordinate, project onto n and axis coordinates
//   coeffs (...); n_coord = sum(coeffs_i * axis_coord^i)
// }
//
// This is an extruded 2D polynomial
//
namespace Foam {
namespace optiMesh {

class polynomial : public optiConstraint
{
  point origin_;

  vector axis_;

  vector n_;

  vector t_;

  scalarField coeffs_;

  // evaluate the polynomial
  scalar evaluate(const scalar& x) const;

  point project(const label&, const point&) const;

public:
  TypeName("polynomial");

  polynomial(const fvMesh&, const dictionary&);

  virtual ~polynomial();

  virtual void constrain(pointField& pf) const;
};

}
}
#endif
