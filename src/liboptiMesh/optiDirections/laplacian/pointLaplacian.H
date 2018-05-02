#ifndef pointLaplacian_H
#define pointLaplacian_H

#include "laplacian.H"

namespace Foam {
namespace optiMesh {

class pointLaplacian : public laplacian
{
public:
  enum weightType {
    simple,
    idw,
    unique // if duplicate points are detected in weighting kernel then give them half the weight
  };

private:
  weightType wType_;

  // tolerance for uniqueness
  scalar tol_;

  scalar calcUniqueness(const label& pointI, const labelList& allPointIds,
    const pointField& points) const;

  virtual scalar weight(const label& pointI, const label& otherI, 
      const labelList& allOtherIds) const;

  virtual void average(const label& pointI, point& pSum, scalar& wSum) const;

public:
  TypeName("pointLaplacian");

  pointLaplacian(const fvMesh&, const dictionary&);

  virtual ~pointLaplacian();

  //virtual void update(); // use parent update function
};

}
}

#endif
