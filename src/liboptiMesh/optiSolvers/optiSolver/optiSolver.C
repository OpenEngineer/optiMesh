#include "optiSolver.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(optiSolver, 0);
  defineRunTimeSelectionTable(optiSolver, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


scalar optiSolver::magSq(const vectorField& a)
{
  scalar sum = 0.0;

  forAll(a, i) {
    sum += a[i]&a[i];
  }

  return sum;
}


scalar optiSolver::magSq(const vectorField& a, const vectorField& b)
{
  scalar sum = 0.0;

  forAll(a, i) {
    sum += a[i] & b[i];
  }

  return sum;
}


optiSolver::optiSolver(const fvMesh& mesh, const dictionary& dict) :
  mesh_(mesh)
{}


optiSolver::~optiSolver()
{}
