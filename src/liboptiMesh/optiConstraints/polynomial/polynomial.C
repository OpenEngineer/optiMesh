#include "polynomial.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(polynomial, 0);
  addToRunTimeSelectionTable(optiConstraint, polynomial, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;

scalar polynomial::evaluate(const scalar& x) const
{
  scalar xx = 1.0;

  scalar sum = 0.0;

  forAll(coeffs_, i) {
    sum += xx*coeffs_[i];
    xx *= x;
  }

  return sum;
}


point polynomial::project(const label& origI, const point& a) const
{
  // dont use the original points

  vector diff = a - origin_;

  // in local coordinate system
  scalar x = diff & axis_;
  scalar y = evaluate(x);
  scalar z = diff & t_;

  point proj = origin_ + x*axis_ + y*n_ + z*t_;

  return proj;
}


polynomial::polynomial(const fvMesh& mesh, const dictionary& dict) : 
  optiConstraint(mesh, dict),
  origin_(dict.lookup("origin")),
  axis_(dict.lookup("axis")),
  n_(dict.lookup("n")),
  t_(axis_ ^ n_),
  coeffs_(dict.lookup("coeffs"))
{
  axis_ /= mag(axis_);
  n_ /= mag(n_);
  t_ /= mag(t_);
}


polynomial::~polynomial()
{}


void polynomial::constrain(pointField& pf) const
{
  forAllConstIter(pointSet, *this, iter) 
  {
    label i = iter.key();

    pf[i] = project(i, pf[i]);
  }
}
