#ifndef none_H
#define none_H

#include "optiSolver.H"

namespace Foam {
namespace optiMesh {

// dont do anything
//
// using this in combination with an objective is equivalent to gradient descent
class none : public optiSolver {
public:
  TypeName("none");

  none(const fvMesh&, const dictionary&);

  virtual ~none();

  virtual void update(vectorField&);
};

}
}
#endif
