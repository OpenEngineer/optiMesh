#ifndef fixed_H
#define fixed_H

#include "optiConstraint.H"

namespace Foam {
namespace optiMesh {

class fixed : public optiConstraint
{
public:
  TypeName("fixed");

  fixed(const fvMesh&, const dictionary&);

  virtual ~fixed();

  virtual void constrain(pointField& pf) const;
};

}
}

#endif
