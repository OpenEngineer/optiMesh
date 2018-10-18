#include "extrudedProfile.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(extrudedProfile, 0);
  addToRunTimeSelectionTable(optiConstraint, extrudedProfile, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


extrudedProfile::extrudedProfile(const fvMesh& mesh, const dictionary& dict) :
  optiConstraint(mesh, dict),
  profile_
  (
    IOField<point>(
      IOobject
      (
        fileName(dict.lookup("profile")),
        mesh.objectRegistry::db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
      )
    )
  )
{
  local2Global_ = labelList(this->size(), -1);

  closest_ = labelList(this->size(), -1);

  label count = 0;
  forAllConstIter(pointSet, *this, iter){
    label pointI = iter.key();

    local2Global_[count] = pointI;

    closest_[count] = naiveClosest(mesh.points()[pointI]);

    count++;
  }
}


extrudedProfile::~extrudedProfile()
{}


point extrudedProfile::projectOntoSegment(point a, point b, point p) const
{
  scalar z = p.z();

  a.z() = 0.0;
  b.z() = 0.0;
  p.z() = 0.0;

  vector n = (b - a);
  scalar dab = Foam::sqrt(n&n);

  // handle degenerate case
  if (dab < SMALL) {
    return a;
  }

  // project onto segment
  n /= dab;

  scalar daproj = (p-a)&n;

  if (daproj < 0.0) {
    return a;
  } else if (daproj > dab) {
    return b;
  } 

  point proj = a + n*daproj;

  proj.z() = z;

  return proj;
}

scalar extrudedProfile::distanceToSegment(point a, point b, point p) const
{
  a.z() = 0.0;
  b.z() = 0.0;
  p.z() = 0.0;

  point proj = projectOntoSegment(a, b, p);

  return Foam::sqrt((p-proj)&(p-proj));
}

label extrudedProfile::naiveClosest(const point& p) const
{
  // For each section calculate the distance
  scalar closestD = GREAT;
  label closestI = -1;

  for(label i = 0; i < profile_.size() - 1; i++) {
    // XXX: ignores the Z-coordinate
    scalar d = distanceToSegment(profile_[i], profile_[i+1], p);

    if (d < closestD) {
      closestD = d;
      closestI = i;
    }
  }

  return closestI;
}

point extrudedProfile::naiveProject(const point& p) const 
{
  label i = naiveClosest(p);

  point proj = projectOntoSegment(profile_[i], profile_[i+1], p);

  return proj;
}

void extrudedProfile::constrain(pointField& pf) const
{
  forAll(local2Global_, i) {
    label pointI = local2Global_[i];

    point p = pf[pointI];
    point pFlat(p.x(), p.y(), 0.0);

    pFlat = naiveProject(pFlat);
    pf[pointI] = point(pFlat.x(), pFlat.y(), p.z());

    // TODO: uses the previously found closest
  }
}
