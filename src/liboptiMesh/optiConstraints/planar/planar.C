#include "planar.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(planar, 0);
  addToRunTimeSelectionTable(optiConstraint, planar, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


point planar::project(const label& origI, const point& a) const
{
  point p = singlePlane_ ? p_  : orig_[origI];

  vector diff = a - p;
  point proj = a - n_*(n_&diff);

  return proj;
}


planar::planar(const fvMesh& mesh, const dictionary& dict) : 
  optiConstraint(mesh, dict),
  n_(dict.lookup("n")),
  singlePlane_(dict.found("p")),
  p_(singlePlane_ ? point(dict.lookup("p")) : point(0.0,0.0,0.0))
{}


planar::~planar()
{}


void planar::constrain(pointField& pf) const
{
  forAllConstIter(pointSet, *this, iter) 
  {
    label i = iter.key();

    pf[i] = project(i, pf[i]);
  }
}
