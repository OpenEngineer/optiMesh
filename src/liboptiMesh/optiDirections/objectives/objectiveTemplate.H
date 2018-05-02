template <class ScalarType>
void objective::penalize(ScalarType& x) const
{
  if (penalizationB_ == 1.0) {
    // dont penalize
    return;
  }

  // first do the mapping
  if (objectiveMapConstant_ != 0.0) {
    x += objectiveMapConstant_;
  }

  if (objectiveMapFactor_ != 1.0) {
    x *= objectiveMapFactor_;
  }

  // quadratic penalization
  ScalarType x2 = x*x;

  x *= penalizationB_;

  x += penalizationA_*x2 + penalizationC_;
}


template <class ScalarType, class PointType>
void objective::faceCentreAndNormal(const label& faceI, const List<PointType>& points,
    const labelList& pointIds, PointType& centre, PointType& normal) const
{
  face f = mesh_.faces()[faceI];

  // first calculate the simple centre
  if (fcType_ == simple || fcType_ == triangles) {
    forAll(f, i) {
      const PointType& p = points[findIndex(pointIds, f[i])];

      centre += p;
    }

    centre /= scalar(f.size());
  }

  List<PointType> ec(f.size());
  List<PointType> pq(f.size());

  // only used in triangles mode
  PointType centre_(0.0,0.0,0.0);
  ScalarType area_(0.0);

  forAll(f, i) {
    const PointType& p = points[findIndex(pointIds, f[i])];
    const PointType& q = points[findIndex(pointIds, f[(i+1)%f.size()])];

    // do slightly different things in the loop depending on the fcType_
    if (fcType_ == simple || fcType_ == gauss) {
      // centre of the edge
      ec[i] = p + q; // defer multiplication by 0.5
      pq[i] = q - p; // i to inext!

      normal += ec[i] ^ pq[i]; // defer multiplication by 0.5 * 0.5
      // (derivation: triArea = 0.5*(ec - centre) ^ pq, but sum over second part is zero)
    } else {
      PointType leg0 = p - centre;
      PointType leg1 = q - centre;

      PointType triNormal = leg0 ^ leg1; // defer division by 2

      normal += triNormal; // defer division by 2

      ScalarType triArea = mag(triNormal); // defer division 2
      PointType triCentre =  p + q + centre; // defer division by 3
      centre_ += triArea*triCentre;
      area_ += triArea; // only used for correct weighting of centre_
    }
  }

  if (fcType_ == simple || fcType_ == gauss) {
    normal *= 0.25; // normal thus has magnitude of area
  } else {
    normal *= 0.5;

    centre_ /= (3.0*area_);

    centre = centre_;
  }


  if (fcType_ == gauss) {
    // more accurate face centre, but just too slow
    // loop a second time to calculate the centre, defer division by area
    forAll(f, i) {
      // calc the outward pointing edge normal, magnitude is length of the edge
      PointType en = (pq[i] ^ normal); // defer division by area * area

      centre += (ec[i] & ec[i]) * en; // defer multiplication by 0.5 * 0.5 * 0.5
    }

    ScalarType areaSq = normal & normal;

    centre *= 0.125 / areaSq;

    // we still need to project the centre onto the face, because we dont have any info along the normal
    centre -= normal * ( (centre - 0.5*ec[0]) & normal) / areaSq; // 0.5*ec[0] is just a reference point
  }

  
  if (debug_) {
    compareToNative(normal, f.normal(mesh_.points()), "faceNormal wrong");

    compareToNative(centre, f.centre(mesh_.points()), "faceCentre wrong");
  }
}


template <class ScalarType, class PointType>
void objective::cellCentreAndVolume(const label& cellI, const List<PointType>& faceCentres,
  const List<PointType>& faceNormals, PointType& centre, ScalarType& volume) const
{
  cell c = mesh_.cells()[cellI];

  forAll(c, i) {
    label faceI = c[i];
      
    scalar alignment = (faceI < mesh_.nInternalFaces() && 
                        mesh_.faceNeighbour()[faceI] == cellI) ? -1.0 : 1.0;

    volume += alignment * (faceCentres[i] & faceNormals[i]); // defer division by 3.0

    centre += alignment * (faceCentres[i] & faceCentres[i]) * faceNormals[i]; // defer multiplication by 0.5 and division by volume
  }

  volume /= 3.0;

  centre *= 0.5 / volume;

  if (debug_) {
    compareToNative(volume, mesh_.cellVolumes()[cellI], "cellVolume wrong");


    /*{
      compareToNative(centre, mesh_.cellCentres()[cellI], "cellCentre wrong");
      vector diff = centre.value() - mesh_.cellCentres()[cellI];
      if (Foam::sqrt(diff&diff) > debugTol_) {
        FatalErrorInFunction << "cellCentre wrong:" << nl << 
          "  native     =" << mesh_.cellCentres()[cellI] << nl <<
          "  chainRule  =" << centre.value() << nl << exit(FatalError);
      }
    }*/
  }
}
