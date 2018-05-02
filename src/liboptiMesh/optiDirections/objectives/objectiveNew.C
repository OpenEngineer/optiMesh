#include "objective.H"


using namespace Foam;
using namespace Foam::optiMesh;


autoPtr<objective> objective::New(const fvMesh& mesh, const dictionary& dict)
{
  const word t(dict.lookup("type"));

  rtstConstructorTable::iterator cstrIter =
    rtstConstructorTablePtr_->find(t);

  if (cstrIter == rtstConstructorTablePtr_->end()) 
  {
    FatalErrorIn("objective::New(const fvMesh&, const dictionary&)")
      << "Unknown objective type "
      << t << nl << nl 
      << "Valid objective types are: " << nl
      << rtstConstructorTablePtr_->sortedToc() << nl
      << exit(FatalError);
  }

  return autoPtr<objective>(cstrIter()(mesh, dict));
}
