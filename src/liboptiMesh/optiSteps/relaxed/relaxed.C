#include "relaxed.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(relaxed, 0);
  addToRunTimeSelectionTable(optiStep, relaxed, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


relaxed::relaxed(const fvMesh& mesh, const dictionary& dict) :
  optiStep(mesh, dict),
  count_(0)
{
  dict.lookup("relaxation") >> relaxation_;

  forAll(relaxation_, i) {
    if (relaxation_[i] < 0.0 || relaxation_[i] > 1.0) {
      FatalErrorInFunction <<
        "invalid relaxation" << exit(FatalError);
    } else if (relaxation_[i] == 0.0) {
      Info << "Warning: relaxation set to 0.0" << endl;
    }
  }
}


relaxed::~relaxed()
{}


scalar relaxed::update(const vectorField&)
{
  scalar r = relaxation_[count_];

  count_++;
  count_ = min(count_, relaxation_.size()-1);

  return r;
}
