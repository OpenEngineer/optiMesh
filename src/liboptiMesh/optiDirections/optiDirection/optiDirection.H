#ifndef optiDirection_H
#define optiDirection_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "vectorField.H"


namespace Foam {
namespace optiMesh {

// propose a direction in which to move the points in order to optimize the mesh
class optiDirection : public vectorField {
protected:
  const fvMesh& mesh_;

public:
  TypeName("optiDirection");

  declareRunTimeSelectionTable
  (
    autoPtr,
    optiDirection,
    rtst,
    (const fvMesh& mesh, const dictionary& dict),
    (mesh, dict)
  );

  optiDirection(const fvMesh&, const dictionary& dict);

  static autoPtr<optiDirection> New(const fvMesh&, const dictionary&);

  virtual ~optiDirection();

  // uses mesh in current position
  virtual void update() = 0; 

  void write() const; // uses a pointMesh as IOobject
};

}
}
#endif

