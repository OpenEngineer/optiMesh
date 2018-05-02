#ifndef orthogonality_H
#define orthogonality_H

#include "objective.H"

namespace Foam {
namespace optiMesh {

class orthogonality : public objective
{
  labelList influencingPoints(const label& faceI) const;

  template <class ScalarType, class PointType>
  void cellCentre(const label& cellI, 
      const List<PointType>& points, const labelList& pointIds,
      const label& faceI, const PointType& faceCentre, const PointType& faceNormal, 
      PointType& centre) const;

  template <class ScalarType, class PointType>
  ScalarType calcFaceOrthogonality(const label& faceI,
      const List<PointType>& points, const labelList& pointIds) const;

public:
  TypeName("orthogonality");

  orthogonality(const fvMesh&, const dictionary&);

  virtual ~orthogonality();

  virtual void update();

  virtual scalar evaluate(const vectorField&);
};

#include "orthogonalityTemplate.H"

}
}

#endif
