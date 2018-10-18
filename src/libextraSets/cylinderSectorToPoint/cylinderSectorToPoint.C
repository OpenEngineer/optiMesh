#include "cylinderSectorToPoint.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
namespace Foam
{
  defineTypeNameAndDebug(cylinderSectorToPoint, 0);
  addToRunTimeSelectionTable(topoSetSource, cylinderSectorToPoint, word);
  addToRunTimeSelectionTable(topoSetSource, cylinderSectorToPoint, istream);
}

const Foam::scalar PI = Foam::constant::mathematical::pi;

using namespace Foam;


topoSetSource::addToUsageTable Foam::cylinderSectorToPoint::usage_
(
  cylinderSectorToPoint::typeName,
  "\n   Usage: cylinderSectorToPoint (originX originY originZ) (axisX axisY axisZ) (refX refY refZ) angle0 angle1\n\n"
  "    Select all points within a cylinder sector, input angles in degrees, ref vector defines zero angle, rh-rule completes the axis system\n\n"
);


void cylinderSectorToPoint::combine(topoSet& set, const bool add) const
{
  const pointField& pf = mesh_.points();

  forAll(pf, pointI) {
    vector d = pf[pointI] - origin_;

    vector r = d - axis_*(axis_ & d); // radially pointing outward vector

    // normalize
    scalar rN = Foam::sqrt(r&r);

    // a point in the middle lies in all sectors
    if (rN > SMALL) {
      r /= rN;

      // in projected space
      scalar x = r & ref_;
      scalar y = r & alt_;

      scalar a = std::atan2(y, x);

      if (a < 0.0) {
        a += PI*2.0;
      }

      while (a < angle0_) {
        a += PI*2.0;
      }

      if (a >= angle0_ && a <= angle1_) {
        addOrDelete(set, pointI, add);
      }
    } else {
      addOrDelete(set, pointI, add);
    }
  }
}

// angles here are already in radians
cylinderSectorToPoint::cylinderSectorToPoint(const polyMesh& mesh, 
    const point& origin,
    const vector& axis,
    const vector& ref,
    const scalar& angle0,
    const scalar& angle1) :
  topoSetSource(mesh),
  origin_(origin),
  axis_(axis),
  ref_(ref),
  angle0_(angle0),
  angle1_(angle1)
{
  axis_ /= Foam::sqrt(axis_&axis_);
  ref_  /= Foam::sqrt(ref_&ref_);
  alt_ = axis_ ^ ref_;
  alt_ /= Foam::sqrt(alt_&alt_);

  ref_ = alt_ ^ axis_;

  while (angle0_ < 0.0) {
    angle0_ += PI*2.0;
  }
  while (angle1_ < angle0_) {
    angle1_ += PI*2.0;
  }
}

cylinderSectorToPoint::cylinderSectorToPoint(const polyMesh& mesh,
    const dictionary& dict)
  :
    cylinderSectorToPoint(mesh, dict.lookup("origin"), dict.lookup("axis"),
        dict.lookup("ref"), readScalar(dict.lookup("angle0"))/180.0*PI,
        readScalar(dict.lookup("angle1"))/180.0*PI)
{}


cylinderSectorToPoint::cylinderSectorToPoint(const polyMesh& mesh,
    const scalar& angle1,
    const scalar& angle0,
    const vector& ref,
    const vector& axis,
    const point& origin) :
  cylinderSectorToPoint(mesh, origin, axis, ref, angle0, angle1)
{}

// is is actually reverse
cylinderSectorToPoint::cylinderSectorToPoint(const polyMesh& mesh, Istream& is) :
  cylinderSectorToPoint(mesh, 
      readScalar(checkIs(is))/180.0*PI, 
      readScalar(checkIs(is))/180.0*PI, 
      vector(checkIs(is)), 
      vector(checkIs(is)), 
      point(checkIs(is)))
{
}


cylinderSectorToPoint::~cylinderSectorToPoint()
{}


void cylinderSectorToPoint::applyToSet
(
  const topoSetSource::setAction action,
  topoSet& set
) const
{
  if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
  {
    Info << "   Adding points in cylinder sector with origin = " << 
      origin_ << ", axis = " << axis_ << endl;

    combine(set, true);
  } else if (action == topoSetSource::DELETE) 
  {
    Info << "   Removing points in cylinder sector with origin = " <<
      origin_ << ", n = " << axis_ << endl;

    combine(set, false);
  }
}
