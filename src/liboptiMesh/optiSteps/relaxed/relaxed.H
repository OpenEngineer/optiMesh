#ifndef relaxed_H
#define relaxed_H

#include "optiStep.H"

namespace Foam {
namespace optiMesh {

// dont do anything
//
// using this in combination with an objective is equivalent to gradient descent
class relaxed : public optiStep {
  label count_;

  // use successive relaxation factors upon each call of step()
  List<scalar> relaxation_;

public:
  TypeName("relaxed");

  relaxed(const fvMesh&, const dictionary&);

  virtual ~relaxed();

  virtual scalar update(const vectorField&);
};

}
}
#endif
