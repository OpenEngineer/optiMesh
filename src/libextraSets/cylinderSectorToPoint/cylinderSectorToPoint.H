#ifndef cylinderSectorToPoint_H
#define cylinderSectorToPoint_H

#include "topoSetSource.H"

namespace Foam
{
class cylinderSectorToPoint : public topoSetSource
{
  static addToUsageTable usage_;

  point origin_;

  vector axis_;

  vector ref_;

  scalar angle0_; // degrees in dicts and input, radians internally

  scalar angle1_; // degrees in dicts and input, radians internally


  vector alt_; // completed by rh-rule from ref and axis

  void combine(topoSet& set, const bool add) const;

public:

  TypeName("cylinderSectorToPoint");

  cylinderSectorToPoint
  (
    const polyMesh& mesh,
    const point& origin,
    const vector& axis,
    const vector& ref,
    const scalar& angle0,
    const scalar& angle1
  );

  cylinderSectorToPoint
  (
    const polyMesh& mesh,
    const dictionary& dict
  );

  cylinderSectorToPoint
  (
    const polyMesh& mesh,
    const scalar& angle1,
    const scalar& angle0,
    const vector& ref,
    const vector& axis,
    const point& origin
  );

  cylinderSectorToPoint
  (
    const polyMesh& mesh,
    Istream&
  );

  virtual ~cylinderSectorToPoint();

  virtual sourceType setType() const
  {
    return POINTSETSOURCE;
  }

  virtual void applyToSet
  (
    const topoSetSource::setAction action,
    topoSet&
  ) const;
};
}
#endif
