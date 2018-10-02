#ifndef halfSpaceToCell_H
#define halfSpaceToCell_H

#include "topoSetSource.H"

namespace Foam
{
class halfSpaceToCell : public topoSetSource
{
  static addToUsageTable usage_;

  point origin_;

  vector n_;

  void combine(topoSet& set, const bool add) const;

public:

  TypeName("halfSpaceToCell");


  halfSpaceToCell
  (
    const polyMesh& mesh,
    const point& origin,
    const vector& n
  );

  halfSpaceToCell
  (
    const polyMesh& mesh,
    const dictionary& dict
  );

  halfSpaceToCell
  (
    const polyMesh& mesh,
    Istream&
  );

  virtual ~halfSpaceToCell();


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
