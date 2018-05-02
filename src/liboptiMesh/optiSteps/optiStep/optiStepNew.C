#include "optiStep.H"


using namespace Foam;
using namespace Foam::optiMesh;


autoPtr<optiStep> optiStep::New(const fvMesh& mesh, const dictionary& dict)
{
  const word t(dict.lookup("type"));

  rtstConstructorTable::iterator cstrIter =
    rtstConstructorTablePtr_->find(t);

  if (cstrIter == rtstConstructorTablePtr_->end()) 
  {
    FatalErrorIn("optiStep::New(const fvMesh&, const dictionary&)")
      << "Unknown optiStep type "
      << t << nl << nl 
      << "Valid optiStep types are: " << nl
      << rtstConstructorTablePtr_->sortedToc() << nl
      << exit(FatalError);
  }

  return autoPtr<optiStep>(cstrIter()(mesh, dict));
}
