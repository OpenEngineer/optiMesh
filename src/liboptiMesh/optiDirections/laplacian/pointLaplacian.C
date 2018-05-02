#include "pointLaplacian.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(pointLaplacian, 0);
  addToRunTimeSelectionTable(optiDirection, pointLaplacian, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


scalar pointLaplacian::calcUniqueness(const label& pointI, const labelList& allPointIds,
    const pointField& points) const
{
  label count = 1;

  point p = points[pointI];

  forAll(allPointIds, i) {
    label otherI = allPointIds[i];

    if (otherI == pointI) {
      continue;
    }

    point otherP = points[otherI];

    if (mag(p-otherP) < tol_) {
      Info << "duplicate kernel point detected" << endl;
      count++;
    }
  }

  scalar w = 1.0/scalar(count);

  if (w < 0.999) { 
    Info << " " << w << endl;
  }
  return 1.0/scalar(count);
}


scalar pointLaplacian::weight(const label& pointI, const label& otherI, const labelList& allOtherIds) const
{
  // always assume simple first
  scalar w(1.0);

  // other weight types
  switch (wType_) {
    case simple:
      break;
    case idw:
      w = calcIDW(mesh_.points()[pointI], mesh_.points()[otherI]);
      break;
    case unique:
      w = calcUniqueness(otherI, allOtherIds, mesh_.points());
      break;
  }

  return w;
}


void pointLaplacian::average(const label& pointI, point& pSum, scalar& wSum) const
{
  const labelList& points = mesh_.pointPoints()[pointI];

  forAll(points, i) {
    label otherI = points[i];

    point p = mesh_.points()[otherI];

    scalar w = weight(pointI, otherI, points);

    pSum += w*p;
    wSum += w;
  }
}


pointLaplacian::pointLaplacian(const fvMesh& mesh, const dictionary& dict) : 
  laplacian(mesh, dict),
  wType_(simple),
  tol_(dict.lookupOrDefault<scalar>("tol", 1e-10))
{
  word wTypeWord(dict.lookupOrDefault<word>("weightType", "simple"));

  if (wTypeWord == "simple") {
    wType_ = simple;
  } else if (wTypeWord == "idw") {
    wType_ = idw;
  } else if (wTypeWord == "unique") {
    wType_ = unique;
  } else {
    FatalErrorInFunction << "Weight type not recognized" << exit(FatalError);
  }
}


pointLaplacian::~pointLaplacian()
{}
