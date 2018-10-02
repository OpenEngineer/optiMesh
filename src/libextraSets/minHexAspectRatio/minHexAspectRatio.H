#ifndef minHexAspectRatio_H
#define minHexAspectRatio_H

#include "topoSetSource.H"

namespace Foam
{
class minHexAspectRatio : public topoSetSource
{
  static addToUsageTable usage_;

  vector greater_;

  vector lesser_;

  scalar minRatio_;

  scalar getDelta(const label& cellI, const scalar& V, vector dir) const;

  void combine(topoSet& set, const bool add) const;

public:

  TypeName("minHexAspectRatio");


  minHexAspectRatio
  (
    const polyMesh& mesh,
    const vector& greater,
    const vector& lesser,
    const scalar& minRatio
  );

  minHexAspectRatio
  (
    const polyMesh& mesh,
    const dictionary& dict
  );

  minHexAspectRatio
  (
    const polyMesh& mesh,
    Istream&
  );

  virtual ~minHexAspectRatio();


  virtual sourceType setType() const
  {
    return CELLSETSOURCE;
  }

  virtual void applyToSet
  (
    const topoSetSource::setAction action,
    topoSet&
  ) const;
};
}

#endif
