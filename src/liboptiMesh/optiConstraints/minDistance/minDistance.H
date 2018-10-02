#ifndef localSmoothin_H
#define localSmoothin_H

#include "optiConstraint.H"

// Synopsis: for each point in p0 and not in p1, find closest initial point in p1
//  then maintain at least a fraction of that distance to that point in p1
//
// system/smoothMeshDict entry:
//
// {
//   type minDistance;
//   set p0;
//   fixedSet p1;
//   fractionOfStart 0.5;
// }
//
namespace Foam {
namespace optiMesh {

class minDistance : public optiConstraint
{
  scalar fractionOfStart_;
  scalar absDist_;

  List<bool> active_; // those in fixedSet are inactive
  List<point> fixedPoints_;
  List<point> startPoints_;
  List<label> local2global_;

public:
  TypeName("minDistance");

  minDistance(const fvMesh&, const dictionary&);

  virtual ~minDistance();

  virtual void constrain(pointField& pf) const;
};

}
}
#endif
