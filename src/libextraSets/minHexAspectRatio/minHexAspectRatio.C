#include "minHexAspectRatio.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
  defineTypeNameAndDebug(minHexAspectRatio, 0);
  addToRunTimeSelectionTable(topoSetSource, minHexAspectRatio, word);
  addToRunTimeSelectionTable(topoSetSource, minHexAspectRatio, istream);
}


using namespace Foam;

topoSetSource::addToUsageTable Foam::minHexAspectRatio::usage_
(
  minHexAspectRatio::typeName,
  "\n   Usage: minHexAspectRatio (originX originY originZ) (nx ny nz)\n\n"
  "    Select all cells on positive side of plane\n\n"
);


scalar minHexAspectRatio::getDelta(const label& cellI, const scalar& V, vector dir) const {
  // Loop the faces,
  const cell& c = mesh_.cells()[cellI];

  label bestFace(-1);
  scalar bestDotProd(0.0);
  scalar bestArea(0.0);

  forAll(c, cellFaceI) {
    label faceI = c[cellFaceI];

    vector faceArea = mesh_.faceAreas()[faceI];
    scalar faceMag = Foam::sqrt(faceArea&faceArea);

    vector n = faceArea/faceMag;

    scalar dotProd = std::abs(n&dir);

    if (dotProd > bestDotProd) {
      bestFace = faceI;
      bestArea = faceMag;
      bestDotProd = dotProd;
    }
  }

  if (bestFace == -1) {
    FatalErrorInFunction << "algo error" << exit(FatalError);
  }

  return V/bestArea;
}

void minHexAspectRatio::combine(topoSet& set, const bool add) const
{
  const scalarField& V = mesh_.cellVolumes();

  forAll(V, cellI)
  {
    scalar dG = getDelta(cellI, V[cellI], greater_);
    scalar dL = getDelta(cellI, V[cellI], lesser_);

    if (scalar(dG/dL) > minRatio_) {
      addOrDelete(set, cellI, add);
    }
  }
}


minHexAspectRatio::minHexAspectRatio(const polyMesh& mesh, const vector& greater, 
    const vector& lesser, const scalar& minRatio) :
  topoSetSource(mesh),
  greater_(greater),
  lesser_(lesser),
  minRatio_(minRatio)
{}


minHexAspectRatio::minHexAspectRatio(const polyMesh& mesh,
    const dictionary& dict)
  :
    minHexAspectRatio(mesh, vector(dict.lookup("greater")), vector(dict.lookup("lesser")), readScalar(dict.lookup("minRatio")))
{}


minHexAspectRatio::minHexAspectRatio
(
  const polyMesh& mesh,
  Istream& is
) :
  topoSetSource(mesh),
  greater_(checkIs(is)),
  lesser_(checkIs(is)),
  minRatio_(readScalar(checkIs(is)))
{}


minHexAspectRatio::~minHexAspectRatio()
{}


void minHexAspectRatio::applyToSet
(
  const topoSetSource::setAction action,
  topoSet& set
) const
{
  if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
  {
    Info << "   Adding cells with " << greater_ << " / " << lesser_ << " aspect ratio larger than " << minRatio_ << endl;

    combine(set, true);
  } else if (action == topoSetSource::DELETE) 
  {
    Info << "   Removing cells with " << greater_ << " / " << lesser_ << " aspect ratio larger than " << minRatio_ << endl;

    combine(set, false);
  }
}
