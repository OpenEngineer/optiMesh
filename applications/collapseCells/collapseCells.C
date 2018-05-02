// collapse cells that have near-zero volume, this means that 4 faces have zero size are to be removed, master face is kept. neighbours that are not marked for collapsing 

#include "fvCFD.H"
#include "polyTopoChange.H"
#include "mergePoints.H"

#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"

labelList unsortedPointMerge(const List<point>& points, scalar tol)
{
  labelList pointMergeMap;
  (void)mergePoints(points, tol, false, pointMergeMap);

  // the pointMergeMap is a sorted one
  labelList comm(points.size(), -1);
  labelList unsorted(points.size(), -1);
  forAll(points, pointI) {
    label commI = pointMergeMap[pointI];

    if (comm[commI] == -1) { // register self
      comm[commI] = pointI;
      unsorted[pointI] = pointI;
    } else {
      unsorted[pointI] = comm[commI];
    }
  }

  // check the result
  label count = 0;
  forAll(unsorted, pointI) {
    if (unsorted[pointI] == -1) {
      FatalErrorInFunction << "unable to fill unsorted pointMergeMap" << exit(FatalError);
    }

    if (unsorted[pointI] != pointI) {
      count++;
    }
  }

  if (count > points.size()/2) {
    FatalErrorInFunction << "unlikely that more than half the points are merged, could be algo error" << exit(FatalError);
  }

  return unsorted;
}


void printCellFaceInfo(label cellI, const fvMesh& mesh, 
    const List<point>& faceCentres, const labelList& faceMergeMap)
{
  cell c = mesh.cells()[cellI];

  Info << nl << "cell " << cellI << " faces: " << endl;

  forAll(c, i) {
    label faceI = c[i];
    label baseI = faceMergeMap[faceI];
    point cf = faceCentres[faceI];
    label ownerI = mesh.faceOwner()[faceI];
    label neighbourI = (faceI < mesh.nInternalFaces()) ? mesh.faceNeighbour()[faceI] : -1;

    Info << "  faceI: " << faceI << ", baseI: " << baseI << nl <<
      "  centre: " << cf << ", ownerI: " << ownerI << ", neighbourI: " << neighbourI << nl <<
      "  V[ownerI]: " << mesh.cellVolumes()[ownerI] << ", V[neighbourI]: " << 
      mesh.cellVolumes()[neighbourI] << endl;
  }
}


void alignFace(face& fNewUnique, const label& ownerI, 
    const label& neighbourI, const fvMesh& mesh)
{
  point c = fNewUnique.centre(mesh.points());
  vector n = fNewUnique.normal(mesh.points());
  n /= mag(n);
  point co = mesh.cellCentres()[ownerI];

  if (neighbourI == -1) {
    vector d = c - co;
    d /= mag(d);

    if ((n&d) < 0.0) {
      fNewUnique.flip();
    }
  } else {
    point cn = mesh.cellCentres()[neighbourI];
    vector d = cn - co;
    d /= mag(d);

    if ((n&d) < 0.0) {
      fNewUnique.flip();
    }
  }
}


// called recursively
void collectConnectedPoints(const label& pointI, const fvMesh& mesh,
    labelList& connectedPointIds, const scalar& tol)
{
  labelList pPoints = mesh.pointPoints()[pointI];

  point p = mesh.points()[pointI];
  forAll(pPoints, i) {
    label otherI = pPoints[i];

    point otherP = mesh.points()[otherI];

    vector diff = p - otherP;

    if (mag(diff) < tol) {
      if (findIndex(connectedPointIds, otherI) == -1) {
        connectedPointIds.append(otherI);
        collectConnectedPoints(otherI, mesh, connectedPointIds, tol);
      }
    }
  }
}


// needed when points accross baffles merge accidentely
labelList mergeCollapsedEdgePoints(const fvMesh& mesh, const scalar& tol)
{
  labelList pointMergeMap(mesh.points().size(), -1);

  forAll(pointMergeMap, pointI) {
    // collect the pointids connect to pointI via its zero length edges
    List<label> connectedPointIds(1, pointI);

    // calls itself recursively
    collectConnectedPoints(pointI, mesh, connectedPointIds, tol);

    sort(connectedPointIds);

    label baseI = connectedPointIds[0];

    pointMergeMap[pointI] = baseI;
  }

  return pointMergeMap;
}


// return pointMergeMap, flag faces for update
labelList collapseEdges(const fvMesh& mesh, const scalar& mergeTol, 
    polyTopoChange& meshChanger, boolList& updateFaceFlags)
{
  labelList pointMergeMap = mergeCollapsedEdgePoints(mesh, mergeTol);

  // remove the points and set any neighbouring faces to update
  pointSet removedPoints(mesh, "removedPoints", IOobject::NO_READ, IOobject::AUTO_WRITE);

  forAll(pointMergeMap, pointI) {
    label baseI = pointMergeMap[pointI];
    if (baseI != pointI) {
      meshChanger.removePoint(pointI, baseI);
      removedPoints.insert(pointI);

      // flag faces for updates
      labelList pointFaces = mesh.pointFaces()[pointI];
      forAll(pointFaces, i) {
        label faceI = pointFaces[i];
        updateFaceFlags[faceI] = true;
      }
    }
  }

  Info << "Writing " << removedPoints.size() << " points to pointSet \"" << 
    removedPoints.name() << "\" " << endl;
  removedPoints.write();

  return pointMergeMap;
}


void collapseCells(const fvMesh& mesh, const scalar& mergeTol,
    polyTopoChange& meshChanger, boolList& updateFaceFlags)
{
  cellSet removedCells(mesh, "removedCells", IOobject::NO_READ, IOobject::AUTO_WRITE);

  // loop the cells of zero volume to make sure those faces are definitely updated
  forAll(mesh.cells(), cellI) {
    if (mesh.cellVolumes()[cellI] < mergeTol) {
      meshChanger.removeCell(cellI, -1);
      removedCells.insert(cellI);

      // mark faces for update
      cell c = mesh.cells()[cellI];
      forAll(c, i) {
        label faceI = c[i];
        updateFaceFlags[faceI] = true;
      }
    }
  }

  Info << "Writing " << removedCells.size() << " cells to cellSet \"" << 
    removedCells.name() << "\" " << endl;
  removedCells.write();
}


List<point> simpleFaceCentres(const fvMesh& mesh)
{
  // merge the faceCentres
  List<point> faceCentres(mesh.faces().size(), point(0.0,0.0,0.0));
  forAll(faceCentres, faceI) {
    face f = mesh.faces()[faceI];
    
    // because of degenerate faces we just take simple average
    forAll(f, i) {
      faceCentres[faceI] += mesh.points()[f[i]];
    }

    faceCentres[faceI] /= scalar(f.size());
  }

  return faceCentres;
}


labelList mergeFaces(const fvMesh& mesh, const scalar& mergeTol)
{
  List<point> faceCentres = simpleFaceCentres(mesh);

  labelList faceMergeMap = unsortedPointMerge(faceCentres, mergeTol);

  return faceMergeMap;
}


// called recursively
void collectConnectedFaces(const label& faceI, const List<point>& faceCentres,
    const label& cellI, const fvMesh& mesh,
    labelList& connectedFaceIds, const scalar& tol)
{
  if (mesh.cellVolumes()[cellI] < tol) {
    point faceP = faceCentres[faceI];

    cell c = mesh.cells()[cellI];

    forAll(c, i) {
      label otherI = c[i];
      point otherP = faceCentres[otherI];

      scalar d = mag(faceP - otherP);

      if (otherI != faceI && d < tol) {
        if (findIndex(connectedFaceIds, otherI) == -1) {
          connectedFaceIds.append(otherI);

          // get other cell
          label ownerI = mesh.faceOwner()[otherI];
          label neighbourI = (otherI < mesh.nInternalFaces()) ? 
            mesh.faceNeighbour()[otherI] : -1;

          if (ownerI == cellI && neighbourI > -1) { // if neighbourI then this is an endpoint (ie. boundary face)
            collectConnectedFaces(otherI, faceCentres, neighbourI, mesh, 
                connectedFaceIds, tol);
          } else {
            collectConnectedFaces(otherI, faceCentres, ownerI, mesh, 
                connectedFaceIds, tol);
          }
        }
      }
    }
  }
}


labelList mergeFaces(const fvMesh& mesh, const List<point>& faceCentres, 
    const scalar& mergeTol)
{
  labelList faceMergeMap(mesh.faces().size(), -1);

  forAll(faceMergeMap, faceI) {
    // collect the pointids connect to pointI via its zero length edges
    List<label> connectedFaceIds(1, faceI);

    // calls itself recursively
    label ownerI = mesh.faceOwner()[faceI];
    label neighbourI = (faceI < mesh.nInternalFaces()) ? 
      mesh.faceNeighbour()[faceI] : -1;

    collectConnectedFaces(faceI, faceCentres, ownerI, 
        mesh, connectedFaceIds, mergeTol);
    if (neighbourI > -1) {
      collectConnectedFaces(faceI, faceCentres, neighbourI, 
          mesh, connectedFaceIds, mergeTol);
    }

    sort(connectedFaceIds);

    label baseI = connectedFaceIds[0];

    faceMergeMap[faceI] = baseI;
  }

  return faceMergeMap;
}


void mergeFaceDetails(const fvMesh& mesh, const labelList& faceMergeMap, 
    const scalar& zeroVolumeTol, labelList& faceMergeOwners,
    labelList& faceMergeNeighbours, labelList& faceMergePatchIds)
{
  forAll(faceMergeMap, faceI) {
    label baseI = faceMergeMap[faceI];
    label ownerI = mesh.faceOwner()[faceI];
    label neighbourI = (faceI < mesh.nInternalFaces()) ? mesh.faceNeighbour()[faceI] : -1;
    label patchI = (neighbourI == -1) ? mesh.boundaryMesh().whichPatch(faceI) : -1;

    // non-zero cellvolume
    if (mesh.cellVolumes()[ownerI] > zeroVolumeTol) {
      if (ownerI != faceMergeOwners[baseI] && ownerI != faceMergeNeighbours[baseI]) {
        if (faceMergeOwners[baseI] == -1) {
          faceMergeOwners[baseI] = ownerI;
        } else if (faceMergeNeighbours[baseI] == -1) {
          faceMergeNeighbours[baseI] = ownerI;
        } else {
          FatalErrorInFunction << "algo error" << exit(FatalError);
        }
      }
    }

    // non-zero cellvolume
    if (neighbourI > -1 && mesh.cellVolumes()[neighbourI] > zeroVolumeTol) {
      if (neighbourI != faceMergeOwners[baseI] && ownerI != faceMergeNeighbours[baseI]) {
        if (faceMergeOwners[baseI] == -1) {
          faceMergeOwners[baseI] = neighbourI;
        } else if (faceMergeNeighbours[baseI] == -1) {
          faceMergeNeighbours[baseI] = neighbourI;
        } else {
          FatalErrorInFunction << "algo error" << exit(FatalError);
        }
      }
    }

    if (patchI > -1) {
      if (faceMergePatchIds[baseI] == -1) {
        faceMergePatchIds[baseI] = patchI;
      } else if (faceMergePatchIds[baseI] != patchI) {
        FatalErrorInFunction << "patch already defined for this face: " << nl <<
         "  faceI: " << faceI << nl << 
         "  patchI: " << faceMergePatchIds[baseI] << nl << 
         "  newPatchI: " << patchI << exit(FatalError);
      }
    }

    // sort
    if (faceMergeOwners[baseI] > faceMergeNeighbours[baseI] && 
      faceMergeNeighbours[baseI] != -1) 
    {
      label tmp = faceMergeOwners[baseI];
      faceMergeOwners[baseI] = faceMergeNeighbours[baseI];
      faceMergeNeighbours[baseI] = tmp;
    }
  }
}


void collapseFaces(const fvMesh& mesh, const scalar& mergeTol,
    const labelList& pointMergeMap, const boolList& updateFaceFlags,
    const scalar& zeroVolumeTol, polyTopoChange& meshChanger)
{
  // merge via zero volume cells instead of the 
  //labelList faceMergeMap = mergeFaces(mesh, zeroVolumeTol);
  labelList faceMergeMap = mergeFaces(mesh, simpleFaceCentres(mesh), zeroVolumeTol);

  faceSet removedFaces(mesh, "removedFaces", IOobject::NO_READ, IOobject::AUTO_WRITE);
  faceSet modifiedFaces(mesh, "modifiedFaces", IOobject::NO_READ, IOobject::AUTO_WRITE);

  // get new neighbours, owners and patchids for the base faces
  labelList faceMergeOwners(faceMergeMap.size(), -1);
  labelList faceMergeNeighbours(faceMergeMap.size(), -1);
  labelList faceMergePatchIds(faceMergeMap.size(), -1);

  mergeFaceDetails(mesh, faceMergeMap, zeroVolumeTol, faceMergeOwners, 
      faceMergeNeighbours, faceMergePatchIds);

  // update the faces (base faces will also be in this list)
  forAll(updateFaceFlags, faceI) {
    if (updateFaceFlags[faceI]) {
      face f = mesh.faces()[faceI];

      face fNew(f);

      forAll(f, i) {
        label pointI = f[i];
        label baseI = pointMergeMap[pointI];

        fNew[i] = baseI;
      }

      // keep the unique points only
      face fNewUnique(1);
      fNewUnique[0] = fNew[0]; // first point is always valid

      // collapse edges of the face
      label count = 0;
      for(label i = 1; i < fNew.size(); i++) {
        if (fNew[i] != fNewUnique[count]) {
          fNewUnique.append(fNew[i]);
          count++;
        }
      }

      // collapse the last edge of the face (connecting last and first vertex)
      if (fNewUnique[count] == fNewUnique[0]) {
        fNewUnique.resize(count); // collapse
      }

      // if fNewUnique is larger than 3, modifyFace, otherwise removeFace
      if (fNewUnique.size() < 3 || faceMergeMap[faceI] != faceI) {
        if (fNewUnique.size() < 3) {
          meshChanger.removeFace(faceI, -1);
        } else {
          meshChanger.removeFace(faceI, faceMergeMap[faceI]);
        }

        removedFaces.insert(faceI);
      } else {
        // faceI == baseI
        label ownerI = faceMergeOwners[faceI];
        label neighbourI = faceMergeNeighbours[faceI];
        label patchI = faceMergePatchIds[faceI];

        if (ownerI == -1) {
          FatalErrorInFunction << "algo error" << exit(FatalError);
        }

        if (neighbourI == -1 && patchI == -1) {
          FatalErrorInFunction << "algo error" << exit(FatalError);
        }

        if (neighbourI != -1 && patchI != -1) {
          FatalErrorInFunction << "algo error" << exit(FatalError);
        }

        // does fNewUnique really point from ownerI to neighbourI
        alignFace(fNewUnique, ownerI, neighbourI, mesh);

        meshChanger.modifyFace(fNewUnique, faceI, ownerI,
            neighbourI, false, patchI, -1, false);
        modifiedFaces.insert(faceI);
      }
    } else if (mesh.faces()[faceI].mag(mesh.points()) < 1e-10) {
      FatalErrorInFunction << "face shouldve been deleted" << exit(FatalError);
    }
  }


  Info << "Writing " << removedFaces.size() << " faces to faceSet \"" << 
    removedFaces.name() << "\" " << endl;
  removedFaces.write();
  Info << "Writing " << modifiedFaces.size() << " faces to faceSet \"" << 
    modifiedFaces.name() << "\" " << endl;
  modifiedFaces.write();
}


int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"


  polyTopoChange::debug = false;
  polyTopoChange meshChanger(mesh);

  boolList updateFaceFlags(mesh.faces().size(), false);


  labelList pointMergeMap = collapseEdges(mesh, 1e-10, meshChanger, updateFaceFlags);

  collapseCells(mesh, 1e-10, meshChanger, updateFaceFlags);

  collapseFaces(mesh, 1e-10, pointMergeMap, updateFaceFlags, 1e-10, meshChanger);

  (void)meshChanger.changeMesh(mesh, false);


  runTime++;


  Info << nl << "Writing collapsed mesh at time " << runTime.timeName() << nl << endl;
  mesh.write();


  Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;

  Info<< "End\n" << endl;


  return 0;
}
