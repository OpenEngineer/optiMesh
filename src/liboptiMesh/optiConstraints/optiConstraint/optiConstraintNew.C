#include "optiConstraint.H"

using namespace Foam;
using namespace Foam::optiMesh;

autoPtr<optiConstraint> optiConstraint::New(const fvMesh& mesh, const dictionary& dict)
{
  const word t(dict.lookup("type"));


  rtstConstructorTable::iterator cstrIter =
    rtstConstructorTablePtr_->find(t);

  if (cstrIter == rtstConstructorTablePtr_->end())
  {
    FatalErrorIn("optiConstraint::New(const fvMesh&, const dictionary&)")
      << "Unknown optiConstraint type "
      << t << nl << nl
      << "Valid optiConstraint types are: " << nl
      << rtstConstructorTablePtr_->sortedToc() << nl
      << exit(FatalError);
  }


  autoPtr<optiConstraint> ptr(cstrIter()(mesh, dict));

  Info << "Selecting optiConstraint type: " << t << endl;

  if (ptr.empty()) {
    FatalErrorInFunction << "ptr cannot be empty here" << exit(FatalError);
  }
  return ptr;
}
