#ifndef faceLaplacian_H
#define faceLaplacian_H

#include "laplacian.H"

namespace Foam {
namespace optiMesh {

class faceLaplacian : public laplacian
{
public:
  enum weightType {
    simple,
    area,
    idw
  };

private:
  weightType wType_;

  virtual scalar weight(const label& pointI, const label& faceI) const;

  virtual void average(const label& pointI, point& pSum, scalar& wSum) const;

public:
  TypeName("faceLaplacian");

  faceLaplacian(const fvMesh&, const dictionary&);

  virtual ~faceLaplacian();

  //virtual void update(); // use parent update function
};

}
}

#endif
