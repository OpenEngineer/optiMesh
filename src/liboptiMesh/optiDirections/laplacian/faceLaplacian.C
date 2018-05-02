#include "faceLaplacian.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(faceLaplacian, 0);
  addToRunTimeSelectionTable(optiDirection, faceLaplacian, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


scalar faceLaplacian::weight(const label& pointI, const label& faceI) const
{
  scalar w(1.0);

  // other weight types
  switch (wType_) {
    case simple:
      break;
    case area:
      {
        face f = mesh_.faces()[faceI];
        w = f.mag(mesh_.points());
      }
      break;
    case idw:
      {
        face f = mesh_.faces()[faceI];
        w = calcIDW(mesh_.points()[pointI], f.centre(mesh_.points()));
      }
      break;
  }

  return w;
}


void faceLaplacian::average(const label& pointI, point& pSum, scalar& wSum) const
{
  const labelList& faces = mesh_.pointFaces()[pointI];

  forAll(faces, i) {
    label faceI = faces[i];

    face f = mesh_.faces()[faceI];

    point p = f.centre(mesh_.points());

    scalar w = weight(pointI, faceI);

    pSum += w*p;
    wSum += w;
  }
}


faceLaplacian::faceLaplacian(const fvMesh& mesh, const dictionary& dict) : 
  laplacian(mesh, dict),
  wType_(simple)
{
  word wTypeWord(dict.lookupOrDefault<word>("weightType", "simple"));

  if (wTypeWord == "simple") {
    wType_ = simple;
  } else if (wTypeWord == "area") {
    wType_ = area;
  } else if (wTypeWord == "idw") {
    wType_ = idw;
  } else {
    FatalErrorInFunction << "Weight type not recognized" << exit(FatalError);
  }
}


faceLaplacian::~faceLaplacian()
{}
