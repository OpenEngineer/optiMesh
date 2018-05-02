#include "fvCFD.H"
#include "wallDist.H"
#include "wallPolyPatch.H"
#include "mathematicalConstants.H"

#include "optiDirection.H"
#include "optiConstraint.H"
#include "optiSolver.H"
#include "optiStep.H"


// TODO: simplify, and apply on a per-patch basis
typedef struct expansionRatio_t {
  scalar r0;
  scalar rInf;
  scalar nLayers; // approximate value
  scalar d0;
  scalar yCut;
  scalar A; // feedback power
} expansionRatio;


expansionRatio initExpansionRatio(const dictionary& dict)
{
  expansionRatio rExp;

  dict.lookup("r0") >> rExp.r0;
  dict.lookup("rInf") >> rExp.rInf;
  dict.lookup("nLayers") >> rExp.nLayers;
  dict.lookup("d0") >> rExp.d0;
  dict.lookup("yCut") >> rExp.yCut;
  dict.lookup("A") >> rExp.A;

  return rExp;
}


void readPositionedPointSets(const fvMesh& mesh, const dictionary& dict, 
    List<autoPtr<optiMesh::optiConstraint>>& positions)
{
  List<dictionary> subDicts(0);
  dict.lookup("constraints") >> subDicts;

  forAll(subDicts, i) {
    positions.append(optiMesh::optiConstraint::New(mesh, subDicts[i]));
  }
}


volScalarField getWallLayers(const fvMesh& mesh, Time& runTime, 
    const expansionRatio& rExp)
{
  volScalarField wl
  (
    IOobject
    (
      "wallLayer",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimless
  );

  // unset value
  wl = -1;
  // ambiguous value is GREAT
  scalar ambiguousVal = rExp.nLayers;

  // loop the mesh polyPatches
  forAll(mesh.boundaryMesh(), patchI) {
    const polyPatch& patch = mesh.boundaryMesh()[patchI];

    Info << "patch: " << patch.name() << endl;
    if (isA<wallPolyPatch>(patch)) {
      Info << "found a wall patch" << endl;
      // give all owner cells 0
      forAll(patch, patchFaceI) {
        label faceI = patchFaceI + patch.start();

        label cellI = mesh.faceOwner()[faceI];

        if (wl[cellI] < -0.5) {
          wl[cellI] = 0.0; // layer 0
        } else { // ambiguous, set by another patch already
          wl[cellI] = ambiguousVal;
        }
      }
    }
  }

  for(label layerI = 1; layerI < rExp.nLayers; layerI++) {
    forAll(mesh.cells(), cellI) {
      // only treat unset cells
      if (wl[cellI] < -0.5) {
        cell c = mesh.cells()[cellI];

        // loop the faces, grabbing the min that is not zero
        forAll(c, cellFaceI) {
          label faceI = c[cellFaceI];

          if (faceI < mesh.faceNeighbour().size()) {
            label neighbourI = mesh.faceNeighbour()[faceI];
            if (neighbourI == cellI) {
              neighbourI = mesh.faceOwner()[faceI];
            }

            if (neighbourI > -1 && 
                std::abs(wl[neighbourI]-(layerI-1.0)) < SMALL)
            {
              if (wl[cellI] > -0.5) {
                wl[cellI] = ambiguousVal;
              } else {
                wl[cellI] = scalar(layerI);
              }
            }
          }
        }
      }
    }
  }

  forAll(wl, i) {
    if (wl[i] < -0.5) {
      wl[i] = ambiguousVal;
    }
  }

  wl.write();

  return wl;
}


scalar wantedYDist(const expansionRatio& rExp, scalar yLayer)
{
  scalar l;

  if (rExp.r0 > 1.0 + SMALL) {
    // distance up to cell
    l = rExp.d0*(Foam::pow(rExp.r0, yLayer)-1.0)/(rExp.r0-1.0);

    // distance in cell
    scalar d = rExp.d0*Foam::pow(rExp.r0, yLayer);

    l += d*0.5;
  } else {
    l = (yLayer+0.5)*rExp.d0;
  }

  return l;
}


point wantedCellCentre(point cc, const expansionRatio& rExp,
    scalar yLayer, scalar yDist, vector yn)
{
  if (yLayer < rExp.nLayers) {
    scalar yWanted = wantedYDist(rExp, yLayer);

    cc -= yn*(yWanted - yDist)*rExp.A;
  }

  return cc;
}


// y is the wall distance field
void smooth(const fvMesh& mesh, 
    pointField& pf, 
    scalar relaxation,
    const expansionRatio& rExp,
    const volScalarField& yLayer, 
    const volScalarField& yDist, 
    const volVectorField& yn)
{
  forAll(pf, pointI) {
    point oldP(pf[pointI]);

    // loop the cells neighbouring the point
    labelList cells(mesh.pointCells()[pointI]);

    point newP(0.0,0.0,0.0);
    scalar vTot(0.0);

    forAll(cells, localCellI) {
      label cellI = cells[localCellI];

      scalar v = mesh.cellVolumes()[cellI];
      point cc = mesh.cellCentres()[cellI];

      // modify cc by wanted 
      cc = wantedCellCentre(cc, rExp, yLayer[cellI], yDist[cellI], yn[cellI]);

      vTot += v;
      newP += v*cc;
    }

    newP /= vTot;

    pf[pointI] = (1.0 - relaxation)*oldP + relaxation*newP;
  }
}


// TODO: boundary cell layer functionality here
void constrain(
    const List<autoPtr<optiMesh::optiConstraint>>& constraints, 
    pointField& pf)
{
  forAll(constraints, i) {
    constraints[i]().constrain(pf);
  }
}


void constrainDir(
    const List<autoPtr<optiMesh::optiConstraint>>& constraints, 
    const pointField& pf, 
    vectorField& optiDir)
{
  pointField pfNew = pf + optiDir;

  constrain(constraints, pfNew);

  optiDir = pfNew - pf;
}


void collapseEdges(const Map<label>& edges, const fvMesh& mesh, pointField& pf)
{
  forAllConstIter(Map<label>, edges, iter)
  {
    label edgeI = iter.key();
    label pointI = iter();

    point pBase = pf[pointI];

    edge e = mesh.edges()[edgeI];

    pf[e.start()] = pBase;
    pf[e.end()] = pBase;
  }
}


bool facesDegenerate(const label& cellI, const label& faceI, const label& faceJ, 
    const fvMesh& mesh, const scalar& degenerateAngle)
{
  
  vector ni = mesh.faces()[faceI].normal(mesh.points());
  ni /= mag(ni);

  if (faceI < mesh.nInternalFaces() && mesh.faceNeighbour()[faceI] == cellI) {
    ni *= -1.0;
  }

  vector nj = mesh.faces()[faceJ].normal(mesh.points());
  nj /= mag(nj);

  if (faceJ < mesh.nInternalFaces() && mesh.faceNeighbour()[faceJ] == cellI) {
    nj *= -1.0;
  }

  if (faceI == faceJ) {
    FatalErrorInFunction << "face ids cannot be equal, algo error" << exit(FatalError);
  }

  bool isDegenerate = (Foam::acos(min(1.0-SMALL,max(-1.0+SMALL,ni&nj))) < degenerateAngle);

  return isDegenerate; 
}


label getFaceFaceCommonEdge(const label& faceI, const label& faceJ, const fvMesh& mesh)
{
  const labelList& faceEdgesI = mesh.faceEdges()[faceI];

  const labelList& faceEdgesJ = mesh.faceEdges()[faceJ];

  forAll(faceEdgesI, i) {
    forAll(faceEdgesJ, j) {
      if (faceEdgesI[i] == faceEdgesJ[j]) {
        return faceEdgesI[i];
      }
    }
  }

  FatalErrorInFunction << "no common edge found" << exit(FatalError);
}


label getPointFaceCommonEdge(const label& pointI, const label& faceI, const label& otherPointI, const fvMesh& mesh)
{
  // loop points of face
  face f = mesh.faces()[faceI];

  label goodOtherPointI = -1;
  forAll(f, i) {
    if (f[i] == pointI) {

      if (f.nextLabel(i) == otherPointI) {
        goodOtherPointI = f.prevLabel(i);
      } else if (f.prevLabel(i) == otherPointI) {
        goodOtherPointI = f.nextLabel(i);
      } else {
        FatalErrorInFunction << "algo error" << exit(FatalError);
      }
    }
  }

  // loop the edges of the face, getting edgeI
  labelList faceEdges = mesh.faceEdges()[faceI];

  forAll(faceEdges, i) {
    label edgeI = faceEdges[i];
    edge faceEdge = mesh.edges()[edgeI];

    if (faceEdge.start() == goodOtherPointI && 
        faceEdge.end() == pointI) 
    {
      return edgeI;
    } else if (faceEdge.start() == pointI &&
               faceEdge.end() == goodOtherPointI)
    {
      return edgeI;
    }
  } 

  FatalErrorInFunction << "algo error" << exit(FatalError);
}


void collapseFaces(const label& cellI, const label& faceI,
    const label& faceJ, const fvMesh& mesh, Map<label>& edges)
{
  // the algorithm here can be a little more expensive, because degenerate cells will be rare
  label commonEdgeI = getFaceFaceCommonEdge(faceI, faceJ, mesh);

  edge commonEdge = mesh.edges()[commonEdgeI];

  label startI = commonEdge.start();

  label endI = commonEdge.end();


  // the edges the have startI and endI in common are the ones to be collapsed
  label edgeStartI = getPointFaceCommonEdge(startI, faceI, endI, mesh);

  label edgeEndI = getPointFaceCommonEdge(endI, faceI, startI, mesh);

  label edgeStartJ = getPointFaceCommonEdge(startI, faceJ, endI, mesh);

  label edgeEndJ = getPointFaceCommonEdge(endI, faceJ, startI, mesh);

  edges.insert(edgeStartI, startI);

  edges.insert(edgeEndI, endI);

  edges.insert(edgeStartJ, startI);

  edges.insert(edgeEndJ, endI);
}


void collapseCells(List<bool>& collapsedCells, Map<label>& edges, const fvMesh& mesh,
    const scalar& degenerateAngle)
{
  // loop all the cells
  forAll(mesh.cells(), cellI) {
    // only collapse once
    if (collapsedCells[cellI]) {
      continue;
    }

    cell c = mesh.cells()[cellI];

    // only for hex cells
    if (c.size() != 6) {
      continue;
    }

    bool cellDegenerate = false;

    // now compare any masked face to any other masked face
    forAll(c, i) {
      forAll(c, j) {
        if (i == j) {
          continue;
        }

        label faceI = c[i];
        label faceJ = c[j];

        if (facesDegenerate(cellI, faceI, faceJ, mesh, degenerateAngle)) {
          // searches for common edge within
          collapseFaces(cellI, faceI, faceJ, mesh, edges);
          cellDegenerate = true;
          break;
        }
      }

      if (cellDegenerate) {
        break;
      }
    }

    if (cellDegenerate) {
      collapsedCells[cellI] = true;
    }
  }
}


void initializeCollapsedEdgesAndCells(const fvMesh& mesh, 
    Map<label>& collapsedEdges, 
    boolList& collapsedCells)
{
  forAll(mesh.edges(), edgeI) {
    edge e = mesh.edges()[edgeI];

    point a = mesh.points()[e.start()];
    point b = mesh.points()[e.end()];

    if (mag(a-b) < 1e-10) {
      collapsedEdges.insert(e.start(), e.end());
    }
  }

  forAll(mesh.cells(), cellI) {
    if (mesh.cellVolumes()[cellI] < 1e-10) {
      collapsedCells[cellI] = true;
    }
  }
}


int main(int argc, char *argv[])
{
  argList::addBoolOption("constant", "dont increment time when saving");
  argList::addOption("dict", "system/optiMeshDict", "alternate dictionary");

#include "setRootCase.H"
#include "createTime.H" 
#include "createMesh.H"

  word dictName = args.optionLookupOrDefault<word>("dict", "system/optiMeshDict");

  // load the dict
  IOdictionary dict 
  (
    IOobject 
    (
      dictName,
      runTime,
      IOobject::MUST_READ,
      IOobject::NO_WRITE,
      false
    )
  );

  // read the global parameters
  label nIters(readLabel(dict.lookup("nIters")));
  label writeInterval(readLabel(dict.lookup("writeInterval")));
  bool constant = args.optionFound("constant");

  // load the objects
  autoPtr<optiMesh::optiDirection> optiDir
  (
    optiMesh::optiDirection::New(mesh, dict.subDict("direction"))
  );

  autoPtr<optiMesh::optiSolver> solver
  (
    optiMesh::optiSolver::New(mesh, dict.subDict("solver"))
  );

  autoPtr<optiMesh::optiStep> step
  (
    optiMesh::optiStep::New(mesh, dict.subDict("step"))
  );

  List<autoPtr<optiMesh::optiConstraint>> constraints(0);
  readPositionedPointSets(mesh, dict, constraints);


  // boundary cell layers: init of data structures
  expansionRatio rExp = initExpansionRatio(dict);
  wallDist wd(mesh);
  volScalarField yLayer = getWallLayers(mesh, runTime, rExp);

  // map edges to pointIds
  bool collapseDegenerateHexCells;
  scalar degenerateAngle(dict.lookupOrDefault<scalar>("degenerateAngle", 0));
  // degrees to radians
  degenerateAngle *= (constant::mathematical::pi/180.0);
  dict.lookup("collapseDegenerateHexCells") >> collapseDegenerateHexCells;
  Map<label> collapsedEdges; // empty construction
  List<bool> collapsedCells(mesh.cells().size(), false);

  // fill with already collapsed edges and cells
  if (collapseDegenerateHexCells) {
    initializeCollapsedEdgesAndCells(mesh, collapsedEdges, collapsedCells);
  }

  bool latestWritten = true;
  for(label iter = 1; iter <= nIters; iter++) {
    Info << "Smooth mesh iter " << iter << endl;

    latestWritten = false;

    pointField pf = mesh.points();

    optiDir().update();

    constrainDir(constraints, pf, optiDir());

    solver().update(optiDir());

    pf += step().update(optiDir()) * optiDir();

    // loop the edges to be collapsed
    collapseEdges(collapsedEdges, mesh, pf);

    constrain(constraints, pf);


    //smooth(mesh, pf, relaxation, rExp, yLayer, wd.y(), wd.n());
    // XXX: add cell layer functionality here

    mesh.movePoints(pf);

    wd.movePoints();

    // update the edges to be collapsed
    if (collapseDegenerateHexCells) {
      collapseCells(collapsedCells, collapsedEdges, mesh, degenerateAngle);
      // but leave the actual collapsing for the next iteration
    }

    if (writeInterval > 0 && (iter%writeInterval) == 0) {
      if (!constant) {
        runTime++;
      }
      mesh.write();
      //optiDir().write();
      latestWritten = true;
    }
  }

  if (!latestWritten) {
    if (!constant) {
      runTime++;
    }
    mesh.write();
    //optiDir().write();
  }
}
