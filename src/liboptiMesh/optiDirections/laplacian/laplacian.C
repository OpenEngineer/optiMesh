#include "laplacian.H"


namespace Foam
{
namespace optiMesh 
{
  defineTypeNameAndDebug(laplacian, 0);
}
}


using namespace Foam;
using namespace optiMesh;


scalar laplacian::calcIDW(const point& p, const point& q) const
{
  vector diff = p - q;
  
  scalar d2 = diff&diff;

  scalar wInv(1.0);

  if (idwPower_ == 2.0) {
    wInv = d2;
  } else if (idwPower_ == 1.0) {
    wInv = Foam::sqrt(d2);
  } else {
    wInv = Foam::pow(d2, 0.5*idwPower_);
  }

  return 1.0 / (wInv + idwTol_);
}


laplacian::laplacian(const fvMesh& mesh, const dictionary& dict) : 
  optiDirection(mesh, dict),
  idwPower_(dict.lookupOrDefault<scalar>("idwPower", 2.0)),
  idwTol_(dict.lookupOrDefault<scalar>("idwTol", SMALL))
{}


laplacian::~laplacian()
{}


void laplacian::update()
{
  // loop all the points
  forAll(*this, pointI) {
    point p = mesh_.points()[pointI];

    point pSum(0.0,0.0,0.0);

    scalar wSum(0.0);

    average(pointI, pSum, wSum);

    point newP = pSum / wSum;

    // remember: save the increment
    this->operator[](pointI) = (newP - p);
  }
}
