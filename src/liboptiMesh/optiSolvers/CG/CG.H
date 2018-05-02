#ifndef CG_H
#define CG_H

#include "optiSolver.H"

namespace Foam {
namespace optiMesh {

// conjugate gradient
class CG : public optiSolver {
  // state
  vectorField rho_;

  scalar prevNumerator_;

public:
  TypeName("CG");

  CG(const fvMesh&, const dictionary&);

  virtual ~CG();

  virtual void update(vectorField&);
};

}
}

#endif
