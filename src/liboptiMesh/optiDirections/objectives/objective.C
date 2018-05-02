#include "objective.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(objective, 0);
  defineRunTimeSelectionTable(objective, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


bool objective::debug_ = false;

scalar objective::debugTol_ = 1e-8;

bool objective::compareToNative(const scalar& opti, const scalar& native, const string& msg)
{
  bool ok = false;

  if (std::abs(opti - native) > debugTol_) {
    FatalErrorIn(msg) << nl <<
      "own:    " << opti << nl <<
      "native: " << native << nl << exit(FatalError);
  } else {
    ok = true;
  }

  return ok;
}

bool objective::compareToNative(const crScalar& opti, const scalar& native, const string& msg)
{
  bool ok = false;

  if (std::abs(opti.value() - native) > debugTol_) {
    FatalErrorIn(msg) << nl <<
      "own:    " << opti.value() << nl <<
      "native: " << native << nl << exit(FatalError);
  } else {
    ok = true;
  }

  return ok;
}

bool objective::compareToNative(const point& opti, const point& native, const string& msg)
{
  vector diff = opti - native;

  bool ok = false;

  if (mag(diff) > debugTol_) {
    FatalErrorIn(msg) << nl <<
      "own   : " << opti << nl <<
      "native: " << native << nl << exit(FatalError);
  } else {
    ok = true;
  }

  return ok;
}

bool objective::compareToNative(const crPoint& opti, const point& native, const string& msg)
{
  vector diff = opti.value() - native;

  bool ok = false;

  if (mag(diff) > debugTol_) {
    FatalErrorIn(msg) << nl <<
      "own   : " << opti.value() << nl <<
      "native: " << native << nl << exit(FatalError);
  } else {
    ok = true;
  }

  return ok;
}


objective::objective(const fvMesh& mesh, const dictionary& dict) :
  objective(mesh, dict, 0.0, 1.0)
{}


objective::objective(const fvMesh& mesh, const dictionary& dict,
    scalar objectiveMapConstant, scalar objectiveMapFactor) :
  optiDirection(mesh, dict),
  fcType_(simple),
  penalizationB_(dict.lookupOrDefault<scalar>("penalization", 1.0)),
  penalizationA_(0.5 - 0.5*penalizationB_),
  penalizationC_(penalizationA_),
  objectiveMapConstant_(objectiveMapConstant),
  objectiveMapFactor_(objectiveMapFactor)
{
  if (penalizationB_ < 1.0) {
    Info << "Warning: a penalization smaller than 1.0 doesnt make much sense" << endl;
  }

  word fcTypeWord(dict.lookupOrDefault<word>("faceCentreType", "simple")); 

  if (fcTypeWord == "simple") {
    fcType_ = simple;
  } else if (fcTypeWord == "triangles") {
    fcType_ = triangles;
  } else if (fcTypeWord == "gauss") {
    fcType_ = gauss;
  } else {
    FatalErrorInFunction << "FaceCentreType not recognized" << exit(FatalError);
  }
}


objective::~objective()
{}


void objective::update()
{
  forAll(*this, pointI) {
    this->operator[](pointI) = vector(0.0,0.0,0.0);
  }
}
