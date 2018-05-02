#ifndef cylindrical_H
#define cylindrical_H

#include "optiConstraint.H"

// Can also be used as a line object
namespace Foam {
namespace optiMesh {

class cylindrical : public optiConstraint
{
  point origin_;

  vector axis_;

  bool singleCylinder_;

  scalar r_;

  scalar radius(const label&) const;
  
  point project(const label&, const point&) const;


public:
  TypeName("cylindrical");

  cylindrical(const fvMesh&, const dictionary&);

  virtual ~cylindrical();

  virtual void constrain(pointField& pf) const;
};

}
}
#endif
