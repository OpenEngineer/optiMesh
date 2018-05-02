#include "argList.H"
#include "fvCFD.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "polyPatch.H"


void removePoints(const fvMesh& mesh,
    polyTopoChange& meshChanger,
    const cellSet& removalCellSet)
{
  forAll(mesh.points(), pointI) {
    const labelList& pointCells = mesh.pointCells()[pointI];

    bool anyNotRemoval = false;
    forAll(pointCells, i) {
      label cellI = pointCells[i];

      if (!removalCellSet.found(cellI)) {
        anyNotRemoval = true;
        break;
      }
    }

    if (!anyNotRemoval) {
      meshChanger.removePoint(pointI, -1);
    }
  }
}


void removeCells(const fvMesh& mesh,
    polyTopoChange& meshChanger,
    const cellSet& removalCellSet)
{
  forAllConstIter(cellSet, removalCellSet, iter) {
    meshChanger.removeCell(iter.key(), -1);
  }
}


void removeFaces(const fvMesh& mesh,
    polyTopoChange& meshChanger,
    const cellSet& removalCellSet,
    const label& exposedPatchId)
{
  forAll(mesh.faces(), faceI) {
    face f = mesh.faces()[faceI];

    label ownerI = mesh.faceOwner()[faceI];
    label neighbourI = (faceI < mesh.nInternalFaces()) ? 
      mesh.faceNeighbour()[faceI] : -1;

    if (removalCellSet.found(ownerI) && 
        (neighbourI == -1 || removalCellSet.found(neighbourI)))
    {
      meshChanger.removeFace(faceI, -1);
    } else if (removalCellSet.found(ownerI) && 
        !removalCellSet.found(neighbourI))
    {
      f.flip();
      meshChanger.modifyFace(f, faceI, neighbourI, -1, false, exposedPatchId, -1, false);
    } else if (!removalCellSet.found(ownerI) &&
        removalCellSet.found(neighbourI)) 
    {
      meshChanger.modifyFace(f, faceI, ownerI, -1, false, exposedPatchId, -1, false);
    }
  }
}


label assertPatch(fvMesh& mesh, const word& patchName)
{
  // create the patch if doesnt already exist
  const polyBoundaryMesh& patches = mesh.boundaryMesh();
  DynamicList<polyPatch*> allPatches;

  label patchId = -1;
  forAll(patches, patchI) {
    const polyPatch& pp = patches[patchI];

    if (pp.name() == patchName) {
      patchId = patchI;
    }

    allPatches.append(pp.clone(patches, allPatches.size(), pp.size(), pp.start()).ptr());
  }

  if (patchId == -1) {
    // create the patch
    patchId = patches.size();

    allPatches.append(polyPatch::New
      (
        word("patch"),
        patchName,
        0,
        mesh.faces().size(),
        patchId,
        patches
      ).ptr());
    allPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(allPatches);
  } else {
    // delete the clones
    forAll(allPatches, i) {
      delete allPatches[i];
    }
  }

  return patchId;
}


int main(int argc, char *argv[])
{
  argList::validArgs.append("cellSet");
  argList::validArgs.append("exposedPatch");

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  word cellSetName(args[1]);

  word exposedName(args[2]);

  label patchId = assertPatch(mesh, exposedName);

  // create the patch name

  polyTopoChange::debug = false;
  polyTopoChange meshChanger(mesh);
  
  cellSet removalCellSet(mesh, cellSetName, IOobject::MUST_READ, IOobject::NO_WRITE);
 
  // remove points that only neighbour removal-cells
  removePoints(mesh, meshChanger, removalCellSet);

  // remove cells of removalCellSet themselves
  removeCells(mesh, meshChanger, removalCellSet);

  // remove faces that only neighbour removal-cells), otherwise modify them to 
  removeFaces(mesh, meshChanger, removalCellSet, patchId);

  (void)meshChanger.changeMesh(mesh, false);

  runTime++;

  Info << nl << "Writing mesh with removed cells at time" << runTime.timeName() << nl << endl;
  mesh.write();

  Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
      << "  ClockTime = " << runTime.elapsedClockTime() << " s"
      << nl << endl;

  Info<< "End\n" << endl;

  return 0;
}
