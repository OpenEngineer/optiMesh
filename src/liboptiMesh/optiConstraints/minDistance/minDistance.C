#include "minDistance.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(minDistance, 0);
  addToRunTimeSelectionTable(optiConstraint, minDistance, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


minDistance::minDistance(const fvMesh& mesh, const dictionary& dict) : 
  optiConstraint(mesh, dict),
  fractionOfStart_(readScalar(dict.lookup("fractionOfStart"))),
  absDist_(dict.lookupOrDefault<scalar>("absDist", 0.0))
{
  // fixed point set
  pointSet fixedSet(mesh, dict.lookup("fixedSet"));

  // init the lists
  active_ = List<bool>(this->size(), true);
  fixedPoints_ = List<point>(this->size());
  startPoints_ = List<point>(this->size());

  local2global_ = List<label>(this->size(), -1);

  bool anyActive(false);

  // fill the lists
  label count(0);
  forAllConstIter(pointSet, *this, iter) 
  {
    label pointI = iter.key();

    local2global_[count] = pointI;
    if (fixedSet[pointI]) {
      active_[count] = false;
    } else {
      anyActive = true;
      point sp = mesh.points()[pointI];
      startPoints_[count] = sp;


      // find the closest fixed point
      label closestI = -1;
      point closestP(0.0, 0.0, 0.0);
      scalar closestD(GREAT);
      forAllConstIter(pointSet, fixedSet, innerIter) {
        label fixedI = innerIter.key();

        point p = mesh.points()[fixedI];

        scalar d = Foam::sqrt((p - sp) & (p - sp)); 
        if (d < closestD) {
          closestI = fixedI;
          closestP = p;
          closestD = d;
        }
      }

      if (closestI == -1) {
        FatalErrorInFunction << "algo error" << exit(FatalError);
      }

      fixedPoints_[count] = closestP;
    }

    count++;
  }

  if (!anyActive) {
    FatalErrorInFunction << "no active points?, check the input" << exit(FatalError);
  }
}


minDistance::~minDistance()
{}


void minDistance::constrain(pointField& pf) const
{
  forAll(local2global_, i) {
    label pointI = local2global_[i];

    if (active_[i]) {
      point sp = startPoints_[i];
      point fp = fixedPoints_[i];
      point tp = pf[pointI];

      scalar oldD = Foam::sqrt((sp - fp)&(sp - fp));
      scalar newD = Foam::sqrt((tp - fp)&(tp - fp));

      if (newD < fractionOfStart_*oldD || newD < absDist_) {
        // constrain using current direction vector

        scalar r = fractionOfStart_*oldD;
        if ((newD < absDist_) && (absDist_ < fractionOfStart_*oldD)) {
          r = absDist_;
        }

        vector n = 0.5*((tp - fp) + (sp - fp)); // use average of both
        n /= Foam::sqrt(n&n); // normalize

        pf[pointI] = fp + n*r;
      }
    }
  }
}
