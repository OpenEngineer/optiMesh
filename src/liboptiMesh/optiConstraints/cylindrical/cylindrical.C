#include "cylindrical.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(cylindrical, 0);
  addToRunTimeSelectionTable(optiConstraint, cylindrical, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


scalar cylindrical::radius(const label& origI) const
{
  point p(orig_[origI]);

  vector diff = p - origin_;
  diff = diff - axis_*(axis_&diff);

  scalar r = Foam::sqrt(diff&diff);

  return r;
}


point cylindrical::project(const label& origI, const point& a) const
{
  scalar r = singleCylinder_ ? r_ : radius(origI);

  vector diff = a - origin_;

  // axial component
  point proj = origin_ + axis_*(axis_&diff);

  // radial component, only to be applied if diff is nonzero
  if (!singleCylinder_ || r_ > SMALL) {
    vector n = diff - axis_*(axis_&diff);
    if (Foam::sqrt(n&n) > SMALL) {
      n /= Foam::sqrt(n&n);
      proj += n*r;
    } else if (r_ > SMALL) {
      Info << "Warning: cylindrical pointSet point in centre of cylindrical with non-zero radius " << r_ << endl;
    }
  } // else use cylindrical as line object

  return proj;
}


cylindrical::cylindrical(const fvMesh& mesh, const dictionary& dict) : 
  optiConstraint(mesh, dict),
  origin_(dict.lookup("origin")),
  axis_(dict.lookup("axis")),
  singleCylinder_(dict.found("r")),
  r_(singleCylinder_ ? readScalar(dict.lookup("r")) : 0.0)
{
  // normalize to be sure
  axis_ /= Foam::sqrt(axis_&axis_);

  if (singleCylinder_ && r_ < 0.0) {
    FatalErrorIn("cylindricalPointSet constructor") <<
      "r must be greater than, or equal, to zero" << exit(FatalError);
  }
}


cylindrical::~cylindrical()
{}


void cylindrical::constrain(pointField& pf) const
{
  forAllConstIter(pointSet, *this, iter) 
  {
    label i = iter.key();

    pf[i] = project(i, pf[i]);
  }
}
