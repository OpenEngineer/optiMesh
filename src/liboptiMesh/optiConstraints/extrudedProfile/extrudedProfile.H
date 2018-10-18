#ifndef extrudedProfile_H
#define extrudedProfile_H

#include "optiConstraint.H"

// system/optiMeshDict entry:
//
// TODO: right now this is a rough 2D approximation: generalize to 3D
// {
//   type extrudedProfile;
//   set p0;
//
//   profile <FILE>; // in OpenFOAM pointField format
// }
//
namespace Foam{
namespace optiMesh{

class extrudedProfile : public optiConstraint
{
  pointField profile_;

  labelList local2Global_;

  labelList closest_; // len of set, closest point on profile, local index in profile_


  point projectOntoSegment(point a, point b, point p) const;

  scalar distanceToSegment(point a, point b, point p) const;

  label naiveClosest(const point& p) const;

  point naiveProject(const point& p) const;


public:
  TypeName("extrudedProfile");

  extrudedProfile(const fvMesh&, const dictionary&);

  virtual ~extrudedProfile();

  virtual void constrain(pointField& pf) const;
};

}
}

#endif
