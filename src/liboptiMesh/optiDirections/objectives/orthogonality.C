#include "orthogonality.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(orthogonality, 0);
  addToRunTimeSelectionTable(optiDirection, orthogonality, rtst);
  addToRunTimeSelectionTable(objective, orthogonality, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


labelList orthogonality::influencingPoints(const label& faceI) const
{
  label ownerI = mesh_.faceOwner()[faceI];
  label neighbourI = mesh_.faceNeighbour()[faceI]; // faceI must be internal face

  labelList points = mesh_.cellPoints()[ownerI];
  points.append(mesh_.cellPoints()[neighbourI]);

  // remove the duplicates
  labelList order;
  Foam::uniqueOrder(points, order);

  labelList uniquePoints(order.size());
  forAll(order, i) {
    uniquePoints[i] = points[order[i]];
  }

  return uniquePoints;
}


orthogonality::orthogonality(const fvMesh& mesh, const dictionary& dict) : 
  objective(mesh, dict, 1.0, 0.5)
{}


orthogonality::~orthogonality()
{}


void orthogonality::update()
{
  // set gradients to zero to start
  objective::update();

  // for debugging purposes
  scalar obj = 0.0;

  label nFaces = mesh_.nInternalFaces();

  // loop the internal faces, calculating the objective per face
  for(label faceI = 0; faceI < nFaces; faceI++) {
    // initialize the points and the map
    labelList pointIds = influencingPoints(faceI);
    label n = pointIds.size();

    // use the static members of crPoint and crScalar to set n
    crScalar::nPoints_ = n;
    crPoint::nPoints_ = n;

    List<crPoint> points(n);
    forAll(points, i) {
      label pointI = pointIds[i];
      points[i] = crPoint(mesh_.points()[pointI], i, n);
    }

    crScalar faceOrth = calcFaceOrthogonality<crScalar, crPoint>(faceI, points, pointIds);
    obj += faceOrth.value();

    // now store the gradients
    forAll(pointIds, i) {
      label pointI = pointIds[i];

      this->operator[](pointI) += faceOrth.gradient(i);
    }

    // reset nPoints_, so that init errors can be caught elsewhere
    crPoint::nPoints_ = 0;
    crScalar::nPoints_ = 0;
  }

  Info << "  orthogonality objective, tot: " << obj << ", avg: " << 
    obj/scalar(nFaces) << endl;
}


scalar orthogonality::evaluate(const vectorField& dir)
{
  scalar obj(0.0);

  label nFaces = mesh_.nInternalFaces();

  // loop the internal faces, calculating the objective per face
  for(label faceI = 0; faceI < nFaces; faceI++) {
    // initialize the points and the map
    labelList pointIds = influencingPoints(faceI);
    label n = pointIds.size();

    List<point> points(n);
    forAll(points, i) {
      label pointI = pointIds[i];
      points[i] = mesh_.points()[pointI] + dir[pointI];
    }

    obj += calcFaceOrthogonality<scalar, point>(faceI, points, pointIds);
  }

  obj /= scalar(nFaces);

  return obj;
}
