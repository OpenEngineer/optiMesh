// reuse some face data
template <class ScalarType, class PointType>
void sphericity::cellAreaAndVolume(const label& cellI, 
    const List<PointType>& points, const labelList& pointIds,
    ScalarType& cellArea, ScalarType& cellVolume) const
{
  cell c = mesh_.cells()[cellI];
  List<PointType> faceCentres(c.size());
  List<PointType> faceNormals(c.size());

  scalar debugFaceAreaSum(0.0);
  forAll(c, i) {
    // init, then calculate
    faceCentres[i] = PointType(0.0,0.0,0.0);
    faceNormals[i] = PointType(0.0,0.0,0.0);

    faceCentreAndNormal<ScalarType, PointType>(c[i], points, pointIds, 
        faceCentres[i], faceNormals[i]);

    cellArea += mag(faceNormals[i]);

    if (debug_) {
      debugFaceAreaSum += mesh_.faces()[c[i]].mag(mesh_.points());
    }
  }

  PointType cellCentre(0.0,0.0,0.0); // dummy
  cellCentreAndVolume<ScalarType, PointType>(cellI, faceCentres, faceNormals, 
      cellCentre, cellVolume);

  if (debug_) {
    compareToNative(cellArea, debugFaceAreaSum, "face area sum wrong");
  }
}


template <class ScalarType, class PointType>
ScalarType sphericity::calcCellSphericity(const label& cellI, 
    const List<PointType>& points, const labelList& pointIds) const
{
    ScalarType cellArea(0.0);
    ScalarType cellVolume(0.0);

    cellAreaAndVolume<ScalarType, PointType>(cellI, points, pointIds, 
        cellArea, cellVolume);

    ScalarType sph = cbrt(cellVolume)/sqrt(cellArea);

    penalize(sph);

    return sph;
}
