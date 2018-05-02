#include "optiSolver.H"


using namespace Foam;
using namespace Foam::optiMesh;


autoPtr<optiSolver> optiSolver::New(const fvMesh& mesh, const dictionary& dict)
{
  const word t(dict.lookup("type"));

  rtstConstructorTable::iterator cstrIter =
    rtstConstructorTablePtr_->find(t);

  if (cstrIter == rtstConstructorTablePtr_->end()) 
  {
    FatalErrorIn("optiSolver::New(const fvMesh&, const dictionary&)")
      << "Unknown optiSolver type "
      << t << nl << nl 
      << "Valid optiSolver types are: " << nl
      << rtstConstructorTablePtr_->sortedToc() << nl
      << exit(FatalError);
  }

  return autoPtr<optiSolver>(cstrIter()(mesh, dict));
}
