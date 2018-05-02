#ifndef objective_H
#define objective_H

#include "optiDirection.H"
#include "chainRule.H"

// the objectives get their own RunTimeSelectionTable, so that they can be used for just the evaluation, the calculation functions should be templated, that scalar and crScalar can be used interchangeably, the evaluation using the scalar and points should then be much faster
namespace Foam {
namespace optiMesh {

// pass the positive gradient as a direction to optimize the mesh
// the objective should be a maximum
class objective : public optiDirection
{
private:
  enum faceCentreType {
    simple,
    triangles,
    gauss
  };

  // only gauss cellCentre supported because it is the easiest to calculate

  faceCentreType fcType_;
  
  // quadratic function penalization mapping
  // once objective has been mapped to between 0 and 1, apply a penalization
  scalar penalizationB_; 

  scalar penalizationA_;

  scalar penalizationC_;


protected:
  // mapping of point/face/cell objective to a value between 0-1
  scalar objectiveMapConstant_;

  scalar objectiveMapFactor_;


  static bool debug_; // compare own calculations to OF calculations
  // (we use gauss rule, OF uses tri/tet decompositions)

  static scalar debugTol_;
      
  // only used for debugging within the template functions
  inline static scalar primitiveValue(const scalar& s)
  { return s; };

  inline static scalar primitiveValue(const crScalar& s)
  { return s.value(); };

  inline static point primitiveValue(const point& p)
  { return p; };

  inline static point primitiveValue(const crPoint& p)
  { return p.value(); };


  static bool compareToNative(const scalar&, const scalar&, const string& msg);

  static bool compareToNative(const crScalar&, const scalar&, const string& msg);

  static bool compareToNative(const point&, const point&, const string& msg);

  static bool compareToNative(const crPoint&, const point&, const string& msg);


  // change obj inplace
  template <class ScalarType>
  void penalize(ScalarType& obj) const;
  
  template <class ScalarType, class PointType>
  void faceCentreAndNormal(const label& faceI, const List<PointType>& points,
      const labelList& pointIds, PointType& centre, PointType& normal) const;

  template <class ScalarType, class PointType>
  void cellCentreAndVolume(const label& cellI, const List<PointType>& faceCentres,
      const List<PointType>& faceNormals, PointType& centre, ScalarType& volume) const;

public:
  TypeName("objective");

  declareRunTimeSelectionTable
  (
    autoPtr,
    objective,
    rtst, // or should we name it 'rtst2' ?
    (const fvMesh& mesh, const dictionary& dict),
    (mesh, dict)
  );

  objective(const fvMesh&, const dictionary&);

  objective(const fvMesh&, const dictionary&, scalar objectiveMapConstant, 
      scalar objectiveMapFactor);

  static autoPtr<objective> New(const fvMesh&, const dictionary&);

  virtual ~objective();

  // sets the optidirection to zero
  virtual void update();

  // evaluate objective at 'currentPosition + d'
  // returned value must be point/face/cell averaged in order to be as mesh-independent as possible
  // gradient stays unaveraged though (absolute value shouldnt matter anyway)
  virtual scalar evaluate(const vectorField& d) = 0;
};

// template functions
#include "objectiveTemplate.H"

}
}

#endif
