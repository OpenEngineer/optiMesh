#ifndef localSmoothin_H
#define localSmoothin_H

#include "optiConstraint.H"

// Synopsis: extra explicit pointLaplacian smoothing after the fact, only wrt to other points in the set
// system/smoothMeshDict entry:
//
// {
//   type localSmoothing;
//   set p0;
//   nIters  10;
//   relaxation 1.0;
// }
//
namespace Foam {
namespace optiMesh {

class localSmoothing : public optiConstraint
{
  label nIters_;
  scalar relaxation_;

  List<label> pointIds_;
  List<bool> active_; // depends on optional localFixedPoints pointSet
  List<List<label>> neighbourIds_;
  List<List<label>> localNeighbourIds_; // indices into pointIds_

public:
  TypeName("localSmoothing");

  localSmoothing(const fvMesh&, const dictionary&);

  virtual ~localSmoothing();

  virtual void constrain(pointField& pf) const;
};

}
}
#endif
