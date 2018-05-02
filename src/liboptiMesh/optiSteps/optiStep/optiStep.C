#include "optiStep.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(optiStep, 0);
  defineRunTimeSelectionTable(optiStep, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


optiStep::optiStep(const fvMesh& mesh, const dictionary& dict) :
  mesh_(mesh)
{}


optiStep::~optiStep()
{}
