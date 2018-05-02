#include "cellLaplacian.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(cellLaplacian, 0);
  addToRunTimeSelectionTable(optiDirection, cellLaplacian, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


scalar cellLaplacian::weight(const label& pointI, const label& cellI) const
{
  scalar w(1.0);

  // other weight types
  switch (wType_) {
    case simple:
      break;
    case volume:
      w = mesh_.cellVolumes()[cellI];
      break;
    case idw:
      w = calcIDW(mesh_.points()[pointI], mesh_.cellCentres()[cellI]);
      break;
  }

  return w;
}


void cellLaplacian::average(const label& pointI, point& pSum, scalar& wSum) const
{
  const labelList& cells = mesh_.pointCells()[pointI];

  forAll(cells, i) {
    label cellI = cells[i];

    point p = mesh_.cellCentres()[cellI];

    scalar w = weight(pointI, cellI);

    pSum += w*p;
    wSum += w;
  }
}


cellLaplacian::cellLaplacian(const fvMesh& mesh, const dictionary& dict) : 
  laplacian(mesh, dict),
  wType_(simple)
{
  word wTypeWord(dict.lookupOrDefault<word>("weightType", "simple"));

  if (wTypeWord == "simple") {
    wType_ = simple;
  } else if (wTypeWord == "volume") {
    wType_ = volume;
  } else if (wTypeWord == "idw") {
    wType_ = idw;
  } else {
    FatalErrorInFunction << "Weight type not recognized" << exit(FatalError);
  }
}


cellLaplacian::~cellLaplacian()
{}
