#ifndef cellLaplacian_H
#define cellLaplacian_H

#include "laplacian.H"

namespace Foam {
namespace optiMesh {

class cellLaplacian : public laplacian
{
public:
  enum weightType {
    simple,
    volume,
    idw
  };

private:
  weightType wType_;

  virtual scalar weight(const label& pointI, const label& faceI) const;

  virtual void average(const label& pointI, point& pSum, scalar& wSum) const;

public:
  TypeName("cellLaplacian");

  cellLaplacian(const fvMesh&, const dictionary&);

  virtual ~cellLaplacian();

  //virtual void update(); // use parent update function
};

}
}

#endif
