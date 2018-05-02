#include "LBFGS.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(LBFGS, 0);
  addToRunTimeSelectionTable(optiSolver, LBFGS, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


label LBFGS::circulate(const label& in) const
{
  label n = rho_.size();

  label out = in%n;

  if (out < 0) {
    out += n;
  }

  return out;
}


void LBFGS::push(const scalar& rho, const vectorField& delta, const vectorField& y)
{
  // ignore the first push request
  if (updateI_ > -1) {
    if (updateI_ < maxSize_) {
      rho_.append(rho);
      delta_.append(delta);
      y_.append(y);
    } else {
      label actualI = circulate(updateI_);

      rho_[actualI] = rho;
      delta_[actualI] = delta;
      y_[actualI] = y;
    }
  }

  updateI_++;
}


vectorField LBFGS::recurse(const label& depthI, const vectorField& x) const
{
  vectorField result(x);

  if (depthI < rho_.size()) {
    label k = circulate(updateI_ - depthI);

    scalar rho = rho_[k];
    scalar deltaX = magSq(delta_[k], x);
    vectorField deltaDelta = rho*delta_[k]*deltaX;

    vectorField deeperRhs = x - rho*y_[k]*deltaX;

    vectorField deeperX = recurse(depthI+1, deeperRhs);

    result = deeperX - rho*delta_[k]*magSq(y_[k], deeperX);
  }

  return result;
}


LBFGS::LBFGS(const fvMesh& mesh, const dictionary& dict) :
  optiSolver(mesh, dict),
  prevU_(mesh.points().size(), vector(0.0,0.0,0.0)),
  prevPoints_(mesh.points()),
  maxSize_(readLabel(dict.lookup("maxSize"))),
  updateI_(-1)
{}


LBFGS::~LBFGS()
{}


void LBFGS::update(vectorField& U)
{
  vectorField delta = mesh_.points() - prevPoints_;

  prevPoints_ = mesh_.points();

  vectorField y = U - prevU_;

  prevU_ = U;

  scalar rho = 1.0/(magSq(delta, y) + SMALL);

  push(rho, delta, y);

  U = recurse(0, U);
}
