#ifndef sphericity_H
#define sphericity_H

#include "objective.H"

namespace Foam {
namespace optiMesh {

class sphericity : public objective
{
  labelList influencingPoints(const label& cellI) const;

  template <class ScalarType, class PointType>
  void cellAreaAndVolume(const label& cellI, 
      const List<PointType>& points, const labelList& pointIds,
      ScalarType& cellArea, ScalarType& cellVolume) const;

  template <class ScalarType, class PointType>
  ScalarType calcCellSphericity(const label& cellI,
      const List<PointType>& points, const labelList& pointIds) const;

public:
  TypeName("sphericity");

  sphericity(const fvMesh&, const dictionary&);

  virtual ~sphericity();

  virtual void update();

  virtual scalar evaluate(const vectorField&);
};

#include "sphericityTemplate.H"

}
}

#endif
