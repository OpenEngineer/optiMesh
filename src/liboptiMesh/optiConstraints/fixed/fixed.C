#include "fixed.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(fixed, 0);
  addToRunTimeSelectionTable(optiConstraint, fixed, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


fixed::fixed(const fvMesh& mesh, const dictionary& dict) : 
  optiConstraint(mesh, dict)
{}


fixed::~fixed()
{}


void fixed::constrain(pointField& pf) const
{
  forAllConstIter(pointSet, *this, iter) {
    label i = iter.key();

    pf[i] = orig_[i];
  }
}
