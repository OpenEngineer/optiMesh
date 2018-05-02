#include "pointMesh.H"
#include "pointFields.H"

#include "optiDirection.H"


namespace Foam
{
namespace optiMesh {
  defineTypeNameAndDebug(optiDirection, 0);
  defineRunTimeSelectionTable(optiDirection, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


optiDirection::optiDirection(const fvMesh& mesh, const dictionary& dict) : 
  vectorField(mesh.points().size(), vector(0.0,0.0,0.0)),
  mesh_(mesh)
{
}


optiDirection::~optiDirection()
{}


void optiDirection::write() const
{
  pointVectorField pf
  (
    IOobject
    (
      "U",
      "0",
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    pointMesh(mesh_),
    dimLength
  );

    //
  // set to zero intially
  forAll(*this, pointI) {
    pf[pointI] = this->operator[](pointI);
  }

  pf.write();
}
