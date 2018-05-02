// reuse some face data
template <class ScalarType, class PointType>
void orthogonality::cellCentre(const label& cellI, 
    const List<PointType>& points, const labelList& pointIds,
    const label& faceI, const PointType& faceCentre, const PointType& faceNormal, 
    PointType& centre) const
{
  ScalarType volume(0.0); // dummy

  cell c = mesh_.cells()[cellI];
  List<PointType> faceCentres(c.size());
  List<PointType> faceNormals(c.size());

  forAll(c, i) {
    if (c[i] == faceI) { // we already calculated before
      faceCentres[i] = faceCentre;
      faceNormals[i] = faceNormal;
    } else {
      // init, then calculate
      faceCentres[i] = PointType(0.0,0.0,0.0);
      faceNormals[i] = PointType(0.0,0.0,0.0);

      faceCentreAndNormal<ScalarType, PointType>(c[i], points, pointIds, 
          faceCentres[i], faceNormals[i]);
    }
  }

  cellCentreAndVolume<ScalarType, PointType>(cellI, faceCentres, faceNormals, 
      centre, volume);
}


template<class ScalarType, class PointType>
ScalarType orthogonality::calcFaceOrthogonality(const label& faceI,
    const List<PointType>& points, const labelList& pointIds) const
{
  PointType faceCentre(0.0,0.0,0.0);
  PointType faceNormal(0.0,0.0,0.0);

  // calculate the face 
  faceCentreAndNormal<ScalarType, PointType>(faceI, points, pointIds, 
      faceCentre, faceNormal); 

  // calculate the chainruled properties of the owner cell
  label ownerI = mesh_.faceOwner()[faceI];
  PointType ownerCentre(0.0,0.0,0.0);
  cellCentre<ScalarType, PointType>(ownerI, points, pointIds, 
      faceI, faceCentre, faceNormal, 
      ownerCentre);

  // calculate the chainruled properties of the neighbour cell
  label neighbourI = mesh_.faceNeighbour()[faceI];
  PointType neighbourCentre(0.0,0.0,0.0);
  cellCentre<ScalarType, PointType>(neighbourI, points, pointIds, 
      faceI, faceCentre, faceNormal, neighbourCentre);

  // normalize the cell to cell vector, and the face normal
  PointType cellToCell = neighbourCentre - ownerCentre;
  cellToCell /= sqrt(cellToCell&cellToCell);

  faceNormal /= sqrt(faceNormal&faceNormal);

  // calculate the actual objective
  ScalarType orth = cellToCell & faceNormal;

  penalize(orth);

  return orth;
}
