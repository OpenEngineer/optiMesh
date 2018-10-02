#include "optiConstraint.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(optiConstraint, 0);
  defineRunTimeSelectionTable(optiConstraint, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


optiConstraint::optiConstraint(const fvMesh& mesh, const dictionary& dict) : 
  pointSet(mesh, word(dict.lookup("set"))),
  mesh_(mesh),
  orig_(mesh.points())
{}


optiConstraint::~optiConstraint()
{}
