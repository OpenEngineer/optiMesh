// Author: Christian Schmitz, christian.schmitz@openengineer.eu

#ifndef quadraticSearch_H
#define quadraticSearch_H

#include "optiStep.H"
#include "autoPtr.H"
#include "objective.H"

namespace Foam {
namespace optiMesh {

// evaluate the objective along a given direction at several positions
// fit a least squares quadratic curve through the evaluated points, and keep moving the presumed optimum
//
// stop upon some convergence criterium
//
// the initialSamples should be close enough to zero
class quadraticSearch : public optiStep {
  autoPtr<objective> obj_; // objective to evaluate at each position, doesnt need to be the same as the optiDirection object (but there probably isnt a good reason not to do so

  List<scalar> initialSamples_; // first sample to take, as a fraction of the shortest edge

  label nMaxIters_;

  scalar tol_; // convergence tolerance, stop when improvement in obj is smaller than this


  scalar shortestEdgeLength(const label& pointI) const;

  scalar scaleDirection(const vectorField& dir) const;

  static void fitQuadratic(const scalarField& x, const scalarField& y,
      scalar& a, scalar& b, scalar &c);

public:
  TypeName("quadraticSearch");

  quadraticSearch(const fvMesh&, const dictionary&);

  virtual ~quadraticSearch();

  virtual scalar update(const vectorField&);
};

}
}

#endif
