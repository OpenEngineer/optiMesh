#include "sphericity.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(sphericity, 0);
  addToRunTimeSelectionTable(optiDirection, sphericity, rtst);
  addToRunTimeSelectionTable(objective, sphericity, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


labelList sphericity::influencingPoints(const label& cellI) const
{
  labelList points = mesh_.cellPoints()[cellI];

  return points;
}


sphericity::sphericity(const fvMesh& mesh, const dictionary& dict) : 
  objective
  (
    mesh, 
    dict, 
    0.0, 
    Foam::sqrt(4.0*Foam::constant::mathematical::pi)/
      Foam::cbrt((4.0/3.0)*Foam::constant::mathematical::pi)
  )
{}


sphericity::~sphericity()
{}


void sphericity::update()
{
  // set gradients to zero to start
  objective::update();

  // for debugging purposes
  scalar obj = 0.0;

  label nCells = mesh_.cells().size();

  // loop the internal faces, calculating the objective per face
  for(label cellI = 0; cellI < nCells; cellI++) {
    // initialize the points and the map
    labelList pointIds = influencingPoints(cellI);
    label n = pointIds.size();

    crScalar::nPoints_ = n;
    crPoint::nPoints_ = n;

    List<crPoint> points(n);
    forAll(points, i) {
      label pointI = pointIds[i];
      points[i] = crPoint(mesh_.points()[pointI], i, n);
    }

    crScalar cellSph = calcCellSphericity<crScalar, crPoint>(cellI, points, pointIds);

    obj += cellSph.value();

    forAll(pointIds, i) {
      label pointI = pointIds[i];

      this->operator[](pointI) += cellSph.gradient(i);
    }
  }

  // reset so that init errors can be caught
  crScalar::nPoints_ = 0;
  crPoint::nPoints_ = 0;

  Info << "  sphericity objective, tot: " << obj << ", avg: " << 
    obj/scalar(nCells) << endl;
}


scalar sphericity::evaluate(const vectorField& dir)
{
  scalar obj = 0.0;

  label nCells = mesh_.cells().size();

  for(label cellI = 0; cellI < nCells; cellI++) {

    // initialize the points and the map
    labelList pointIds = influencingPoints(cellI);
    label n = pointIds.size();

    List<point> points(n);
    forAll(points, i) {
      label pointI = pointIds[i];
      points[i] = mesh_.points()[pointI] + dir[pointI];
    }

    scalar cellSph = calcCellSphericity<scalar, point>(cellI, points, pointIds);

    obj += cellSph;
  }
  
  obj /= scalar(nCells);

  return obj;
}
