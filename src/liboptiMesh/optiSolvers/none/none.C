#include "none.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(none, 0);
  addToRunTimeSelectionTable(optiSolver, none, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


none::none(const fvMesh& mesh, const dictionary& dict) :
  optiSolver(mesh, dict)
{}


none::~none()
{}


void none::update(vectorField&)
{
  // dont do anything
}
