#ifndef laplacian_H
#define laplacian_H

#include "optiDirection.H"

namespace Foam {
namespace optiMesh {

// abstract class
class laplacian : public optiDirection
{
  scalar idwPower_; // optional 

  scalar idwTol_; // small number so we dont get division by zero

  virtual void average(const label& pointI, point& pSum, scalar& wSum) const = 0;

protected:
  scalar calcIDW(const point& p, const point& q) const;

public:
  TypeName("laplacian");

  laplacian(const fvMesh&, const dictionary&);

  virtual ~laplacian();

  virtual void update();
};

}
}

#endif
