#ifndef LBFGS_H
#define LBFGS_H

#include "optiSolver.H"

namespace Foam
{
namespace optiMesh 
{

class LBFGS : public optiSolver
{
  vectorField prevU_; // store previous output
  pointField  prevPoints_; // store previous positions

  label maxSize_;
  label updateI_;

  // keep appending to these lists untill the maxSize_ is reached
  List<scalar> rho_;
  List<vectorField> delta_;
  List<vectorField> y_;

  // input index can be negative or larger than maxSize,
  //   output index is always in valid range
  label circulate(const label&) const;

  void push(const scalar& rho, const vectorField& delta, const vectorField& y);
  
  // multiply recursively by LBFGS approximations of inverse Hessian
  vectorField recurse(const label& depthI, const vectorField& x) const;
 
public:
  TypeName("LBFGS");

  LBFGS(const fvMesh&, const dictionary&);

  virtual ~LBFGS();

  virtual void update(vectorField&);
};

}
}
#endif
