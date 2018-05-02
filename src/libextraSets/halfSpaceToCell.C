#include "halfSpaceToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
  defineTypeNameAndDebug(halfSpaceToCell, 0);
  addToRunTimeSelectionTable(topoSetSource, halfSpaceToCell, word);
  addToRunTimeSelectionTable(topoSetSource, halfSpaceToCell, istream);
}


using namespace Foam;

topoSetSource::addToUsageTable Foam::halfSpaceToCell::usage_
(
  halfSpaceToCell::typeName,
  "\n   Usage: halfSpaceToCell (originX originY originZ) (nx ny nz)\n\n"
  "    Select all cells on positive side of plane\n\n"
);


void halfSpaceToCell::combine(topoSet& set, const bool add) const
{
  const pointField& ctrs = mesh_.cellCentres();

  forAll(ctrs, cellI)
  {
    vector d = ctrs[cellI] - origin_;

    // project onto normal
    scalar proj = d & n_;

    if (proj > 0.0) {
      addOrDelete(set, cellI, add);
    }
  }
}


halfSpaceToCell::halfSpaceToCell(const polyMesh& mesh, const point& origin, 
    const vector& n) :
  topoSetSource(mesh),
  origin_(origin),
  n_(n)
{}


halfSpaceToCell::halfSpaceToCell(const polyMesh& mesh,
    const dictionary& dict)
  :
    halfSpaceToCell(mesh, dict.lookup("origin"), dict.lookup("n"))
{}


halfSpaceToCell::halfSpaceToCell
(
  const polyMesh& mesh,
  Istream& is
) :
  topoSetSource(mesh),
  origin_(checkIs(is)),
  n_(checkIs(is))
{}


halfSpaceToCell::~halfSpaceToCell()
{}


void halfSpaceToCell::applyToSet
(
  const topoSetSource::setAction action,
  topoSet& set
) const
{
  if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
  {
    Info << "   Adding cells with centre on positive side of plane with origin = " << 
      origin_ << ", n = " << n_ << endl;

    combine(set, true);
  } else if (action == topoSetSource::DELETE) 
  {
    Info << "   Removing cells with centre on positive side of plane with origin = " <<
      origin_ << ", n = " << n_ << endl;

    combine(set, false);
  }
}
