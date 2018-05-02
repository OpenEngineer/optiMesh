#include "chainRule.H"

using namespace Foam;
using namespace Foam::optiMesh;

int main(int argc, char* argv[]) 
{
  Info << "Testing primitive chain-rules calculations..." << endl;
  // test every operation that can be performed on two points of a line
  //
  // 0. setup points
  crPoint p0(1,0,0,0,3);
  crPoint p1(0,1,0,1,3);
  crPoint p2(0,0,1,2,3);

  // 1. check length
  crScalar d = sqrt((p1 - p0) & (p1 - p0));
  Info << "1. line length " << endl;
  Info << " length of line: " << d.value() << endl;
  Info << " gradient wrt 0: " << d.gradient(0) << endl;
  Info << " gradient wrt 1: " << d.gradient(1) << endl;

  // 2. check length to sum of points
  crScalar dO = sqrt((p1 + p0) & (p1 + p0));
  Info << "2. distance of sum to origin " << endl;
  Info << " distance      : " << dO.value() << endl;
  Info << " gradient wrt 0: " << dO.gradient(0) << endl;
  Info << " gradient wrt 1: " << dO.gradient(1) << endl;

  // 3. normalization of length
  crPoint nO = (p1 + p0) / dO;
  crScalar dnO = sqrt(nO&nO);
  Info << "3. normal from origin to sum " << endl;
  Info << " n             : " << nO.value() << endl;
  Info << " gradient wrt 0: " << nO.gradient(0) << endl;
  Info << " gradient wrt 1: " << nO.gradient(1) << endl;

  Info << "4. length of normal " << endl;
  Info << " dnO           : " << dnO.value() << endl;
  Info << " gradient wrt 0: " << dnO.gradient(0) << endl;
  Info << " gradient wrt 1: " << dnO.gradient(1) << endl;

  crPoint p2Scaled = p2*dnO*dnO*dnO/dnO/dnO*(dnO+dnO)/(dnO-dnO+dnO+dnO);
  Info << "5. scaled point " << endl;
  Info << " value         : " << p2Scaled.value() << endl;
  Info << " gradient wrt 0: " << p2Scaled.gradient(0) << endl;
  Info << " gradient wrt 1: " << p2Scaled.gradient(1) << endl;
  Info << " gradient wrt 2: " << p2Scaled.gradient(2) << endl;

  return 0;
}
