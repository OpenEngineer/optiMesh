#include "optiDirection.H"


using namespace Foam;
using namespace Foam::optiMesh;


autoPtr<optiDirection> optiDirection::New(const fvMesh& mesh, const dictionary& dict)
{
  const word t(dict.lookup("type"));

  rtstConstructorTable::iterator cstrIter =
    rtstConstructorTablePtr_->find(t);

  if (cstrIter == rtstConstructorTablePtr_->end()) 
  {
    FatalErrorIn("optiDirection::New(const fvMesh&, const dictionary&)")
      << "Unknown optiDirection type "
      << t << nl << nl 
      << "Valid optiDirection types are: " << nl
      << rtstConstructorTablePtr_->sortedToc() << nl
      << exit(FatalError);
  }

  return autoPtr<optiDirection>(cstrIter()(mesh, dict));
}
