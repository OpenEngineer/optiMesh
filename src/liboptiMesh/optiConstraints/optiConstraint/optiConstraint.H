#ifndef optiConstraint_H
#define optiConstraint_H

#include "pointSet.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "pointFields.H"

namespace Foam {
namespace optiMesh {

class optiConstraint : public pointSet 
{
protected:
  pointField orig_;

public:
  TypeName("optiConstraint");

  declareRunTimeSelectionTable
  (
    autoPtr,
    optiConstraint,
    rtst,
    (const fvMesh& mesh, const dictionary& dict),
    (mesh, dict)
  );

  optiConstraint(const fvMesh&, const dictionary&);

  static autoPtr<optiConstraint> New(const fvMesh&, const dictionary&);
  
  virtual ~optiConstraint();

  virtual void constrain(pointField& pf) const = 0;
};

}
}

#endif
