#include "CG.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(CG, 0);
  addToRunTimeSelectionTable(optiSolver, CG, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


CG::CG(const fvMesh& mesh, const dictionary& dict) :
  optiSolver(mesh, dict),
  rho_(mesh.points().size()),
  prevNumerator_(1.0)
{
  forAll(rho_, i) {
    rho_[i] = vector(0.0,0.0,0.0);
  }
}


CG::~CG()
{}


// U is a gradient in the direction of an objective that needs to be maximized
void CG::update(vectorField& U)
{
  scalar numerator = magSq(U);

  scalar denominator = prevNumerator_;

  scalar cgFactor = numerator / denominator;

  prevNumerator_ = numerator;

  forAll(rho_, i) {
    rho_[i] = U[i] + rho_[i]*cgFactor;
    
    // store
    U[i] = rho_[i];
  }
}
