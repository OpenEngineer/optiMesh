#include "localSmoothing.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace optiMesh
{
  defineTypeNameAndDebug(localSmoothing, 0);
  addToRunTimeSelectionTable(optiConstraint, localSmoothing, rtst);
}
}


using namespace Foam;
using namespace Foam::optiMesh;


localSmoothing::localSmoothing(const fvMesh& mesh, const dictionary& dict) : 
  optiConstraint(mesh, dict),
  nIters_(readLabel(dict.lookup("nIters"))),
  relaxation_(readScalar(dict.lookup("relaxation")))
{
  // list the pointIds and the pointNeighbours
  pointIds_ = List<label>(this->size());
  label count(0);
  forAllConstIter(pointSet, *this, iter) 
  {
    label i = iter.key();

    pointIds_[count] = i;

    count++;
  }

  // temporary global2local indexing (for easier algorithms)
  List<label> global2local(mesh_.points().size(), -1);
  forAll(pointIds_, i) {
    global2local[pointIds_[i]] = i;
  }

  // optionally set some points to inactive
  active_ = List<bool>(this->size(), true);
  if (dict.found("localFixedPoints")) {
    pointSet lfp(mesh, word(dict.lookup("localFixedPoints")));

    forAllConstIter(pointSet, lfp, iter) {
      label i = iter.key();

      if (this->operator[](i)) {
        active_[global2local[i]] = false;
      }
    }
  }

  // list the neighbourIds of the pointIds, that lie in the pointSet
  neighbourIds_ = List<List<label>>(this->size());
  forAll(pointIds_, i) {
    label pointI = pointIds_[i];

    // gradually append
    List<label> neighbours(0);

    const labelList& allOtherPoints = mesh_.pointPoints()[pointI];
    forAll(allOtherPoints, j) {
      if ((pointI != allOtherPoints[j]) && this->operator[](allOtherPoints[j])) {
        neighbours.append(allOtherPoints[j]);
      }
    }

    if (neighbours.size() == 0) {
      FatalErrorInFunction << "no neighbours detected, algo error" << exit(FatalError);
    }

    neighbourIds_[i] = neighbours;
  }

  // local neighbour indexing 
  localNeighbourIds_ = List<List<label>>(this->size());
  forAll(pointIds_, i) {
    List<label> localNeighbours(neighbourIds_[i].size());
    forAll(localNeighbours, j) {
      label localI = global2local[neighbourIds_[i][j]];
      if (localI == -1) {
        FatalErrorInFunction << "neighbour not in global2local, algo error" << exit(FatalError);
      }

      localNeighbours[j] = localI;
    }

    localNeighbourIds_[i] = localNeighbours;
  }
}


localSmoothing::~localSmoothing()
{}


void localSmoothing::constrain(pointField& pf) const
{
  // init the positions
  List<point> positions(pointIds_.size());
  List<point> newPositions(pointIds_.size());
  forAll(pointIds_, i) {
    label pointI = pointIds_[i];

    positions[i] = pf[pointI];
  }

  for(label iter = 0; iter < nIters_; iter++) {
    forAll(pointIds_, i) {

      scalar wSum(0.0);
      point  pSum(0.0,0.0,0.0);
      forAll(localNeighbourIds_[i], j) {
        scalar w = 1.0;
        point p = positions[localNeighbourIds_[i][j]];
        pSum += p*w;
        wSum += w;
      }
      pSum /= wSum;

      newPositions[i] = relaxation_*pSum + (1.0-relaxation_)*positions[i];
    }

    // copy into positions (only if active though)
    forAll(pointIds_, i) {
      if (active_[i]) {
        positions[i] = newPositions[i];
      }
    }
  }

  // write back to pointField
  forAll(pointIds_, i) {
    label pointI = pointIds_[i];

    pf[pointI] = positions[i];
  }
}
