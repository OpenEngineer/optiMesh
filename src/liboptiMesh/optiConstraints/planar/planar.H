#ifndef planar_H
#define planar_H

#include "optiConstraint.H"

// system/smoothMeshDict entry:
//
// {
//   type planar;
//   set p0;
//   n (0 0 1);
//   p (0 0 0); // optional, otherwise planar to original point
// }
//
namespace Foam {
namespace optiMesh {

class planar : public optiConstraint
{
  vector n_;

  bool singlePlane_; // defaults to false

  point p_; // single planar point

  point project(const label&, const point&) const;

public:
  TypeName("planar");

  planar(const fvMesh&, const dictionary&);

  virtual ~planar();

  virtual void constrain(pointField& pf) const;
};

}
}
#endif
