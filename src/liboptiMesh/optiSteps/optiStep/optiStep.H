#ifndef optiStep_H
#define optiStep_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "vectorField.H"


namespace Foam {
namespace optiMesh {

// modify the optiDirection in order to speed up convergence
class optiStep {
protected:
  const fvMesh& mesh_;

public:
  TypeName("optiStep");

  declareRunTimeSelectionTable
  (
    autoPtr,
    optiStep,
    rtst,
    (const fvMesh& mesh, const dictionary& dict),
    (mesh, dict)
  );

  optiStep(const fvMesh&, const dictionary& dict);

  static autoPtr<optiStep> New(const fvMesh&, const dictionary&);

  virtual ~optiStep();

  virtual scalar update(const vectorField&) = 0;
};

}
}

#endif
