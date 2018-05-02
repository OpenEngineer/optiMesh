// Author: Christian Schmitz, christian.schmitz@openengineer.eu

#include "quadraticSearch.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(quadraticSearch, 0);
  addToRunTimeSelectionTable(optiStep, quadraticSearch, rtst);
}
}


using namespace Foam;
using namespace optiMesh;


// length of shortest edge departing or arriving at pointI
// needs to be recalculated every time because the mesh moves
scalar quadraticSearch::shortestEdgeLength(const label& pointI) const
{
  point p = mesh_.points()[pointI];

  const labelList& otherPointIds = mesh_.pointPoints()[pointI];

  scalar l = GREAT;

  forAll(otherPointIds, i) {
    label otherI = otherPointIds[i];

    point otherP = mesh_.points()[otherI];

    scalar lTest = Foam::mag(p - otherP);

    l = min(lTest, l);
  }

  return l;
}


scalar quadraticSearch::scaleDirection(const vectorField& dir) const
{
  scalar sMin = GREAT;

  forAll(mesh_.points(), pointI) {
    scalar l = shortestEdgeLength(pointI);
    scalar d = Foam::mag(dir[pointI]);

    scalar s = l/(d + SMALL); // +SMALL to avoid division by 0.0

    sMin = min(s, sMin);
  }

  return sMin; // if the input direction is rescaled by this amount, than no point will move further than its closest neighbour (unless two points move towards eachother ofcourse)
}


void quadraticSearch::fitQuadratic(const scalarField& x, 
    const scalarField& y, scalar& a, scalar& b, scalar& c)
{
  // least squares fit
  //
  // The system to solve is: (X^T X) Beta = X^T y
  //
  // we can use a tensor for the LHS matrix
  scalarField Xa(x*x);
  const scalarField& Xb(x);
  scalarField Xc(x.size(), 1.0);

  vector rhs
  (
    sumProd(Xa, y),
    sumProd(Xb, y),
    sumProd(Xc, y)
  );

  tensor lhs
  (
    sumProd(Xa, Xa), sumProd(Xa, Xb), sumProd(Xa, Xc),
    sumProd(Xb, Xa), sumProd(Xb, Xb), sumProd(Xb, Xc),
    sumProd(Xc, Xa), sumProd(Xc, Xb), sumProd(Xc, Xc)
  );

  vector beta = lhs.inv() & rhs;

  a = beta[0];
  b = beta[1];
  c = beta[2];
}


quadraticSearch::quadraticSearch(const fvMesh& mesh, const dictionary& dict) :
  optiStep(mesh, dict),
  obj_(objective::New(mesh, dict.subDict("objective"))),
  initialSamples_(dict.lookup("initialSamples")),
  nMaxIters_(readLabel(dict.lookup("nMaxIters"))),
  tol_(readScalar(dict.lookup("tol")))
{
  if (initialSamples_.size() == 0) {
    initialSamples_.append(0.1);
    initialSamples_.append(0.2);
    Info << "Warning: 0 initial samples specified, using " <<
     initialSamples_[0] << " and " << initialSamples_[1] << endl;
  } else if (initialSamples_.size() == 1) {
    initialSamples_.append(0.5*initialSamples_[0]);
    Info << "Warning: only 1 initial sample specified, using the specified value (" << 
      initialSamples_[0] << " and half the specified value " << initialSamples_[1] << endl;
  }


  bool zeroSpecified = false;
  forAll(initialSamples_, i) {
    if (initialSamples_[i] == 0.0) {
      zeroSpecified = true;
    }
  }

  if (!zeroSpecified) {
    initialSamples_.append(0.0);
  }

  Foam::sort(initialSamples_);
}


quadraticSearch::~quadraticSearch()
{}


scalar quadraticSearch::update(const vectorField& dir)
{
  // the input vector field is just a direction, it contains no step size info
  // this means that it can be very large or very small. we need to rescale it first to lie in usefull range
  scalar s = scaleDirection(dir);

  // perform the initial samples
  scalarField x(initialSamples_.size());
  scalarField y(initialSamples_.size());

  forAll(initialSamples_, i) {
    x[i] = initialSamples_[i]*s;
    y[i] = obj_().evaluate(x[i]*dir);
  }

  // best up till now
  label iPrev = Foam::findMax(y);
  scalar xPrev = x[iPrev];
  scalar yPrev = y[iPrev];

  for (label iterI = 0; iterI < nMaxIters_; iterI++) {
    Info <<  "  " << iterI << ", step size: "  << xPrev << ", obj: " << yPrev << endl;
    // do a quadratic fit of y=a*x^2 + b*x + c, using least squares
    scalar a, b, c;
    fitQuadratic(x, y, a, b, c);

    if (a == 0.0) {
      Info << "Warning in quadraticSearch: zero second derivative" << endl;

      // quit immediately
      // TODO: sample somewhere else?
      break;
    } else if (a > 0.0) {
      Info << "Warning in quadraticSearch: saddle point detected" << endl;

      // quit immediately
      // TODO: sample somewhere else?
      break;
    } else {
      scalar xNext = -b/(2.0*a);
      scalar yNext = obj_().evaluate(xNext*dir);

      bool converged = (yNext  < yPrev + tol_); // no, or small, improvement

      if (yNext > yPrev) { // only store improvements in (xPrev, yPrev)
        xPrev = xNext;
        yPrev = yNext;
      }

      if (converged) {
        break;
      }

      x.append(xNext);
      y.append(yNext);
    }
  }

  return xPrev;
}
