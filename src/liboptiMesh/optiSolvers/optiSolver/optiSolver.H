#ifndef optiSolver_H
#define optiSolver_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "vectorField.H"


namespace Foam {
namespace optiMesh {

// modify the optiDirection in order to speed up convergence
class optiSolver {
protected:
  const fvMesh& mesh_;

  static scalar magSq(const vectorField&);

  static scalar magSq(const vectorField&, const vectorField&);

public:
  TypeName("optiSolver");

  declareRunTimeSelectionTable
  (
    autoPtr,
    optiSolver,
    rtst,
    (const fvMesh& mesh, const dictionary& dict),
    (mesh, dict)
  );

  optiSolver(const fvMesh&, const dictionary& dict);

  static autoPtr<optiSolver> New(const fvMesh&, const dictionary&);

  virtual ~optiSolver();

  virtual void update(vectorField&) = 0;
};

}
}

#endif
