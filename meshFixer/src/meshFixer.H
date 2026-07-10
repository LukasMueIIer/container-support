#ifndef meshFixer_H
#define meshFixer_H

#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "IOdictionary.H"
#include "boolList.H"

#include "motionSmoother.H"
#include "pointMesh.H"
#include "indirectPrimitivePatch.H"

#include "Map.H"
#include "DynamicList.H"
#include "Random.H"

namespace Foam
{

class meshFixer
{
    fvMesh& mesh_;
    Time& runTime_;

    // Workflow control
    bool haveBadFaces_;
    bool discreteStep_;
    bool laplacianStep_;
    bool wiggleStep_;

    // Dictionary parameters
    label totalIterations_;
    bool writeIntermediate_;
    bool writeDisplacement_;
    scalar localRelaxation_;

    label nSmoothLayers_;
    label laplacianIters_;
    scalar localLaplacianRelaxation_;

    label movementRelaxIters_;
    scalar movementRelaxation_;

    scalar wVol_;
    scalar wLap_;
    scalar wExp_;
    scalar maxExpRatio_;
    scalar maxARShield_;

    scalar qGood_; 
    scalar qBad_;  
    scalar sGood_; 
    scalar sBad_;  
    scalar diffConstant_;
    scalar orthoSpring_;
    scalar skewSpring_;

    Random rndGen;

    dictionary qualityControlDict_;
    dictionary relativeQualityControlDict_;
    dictionary wiggleQualityControlDict_;

    // Point mask for locked boundaries
    boolList isBoundaryPoint_;
    boolList isBadPoint_;

    // Private methods
    void identifyBoundaries();
    void calcLocalCorrections(pointVectorField& localDisp);
    void solveGlobalDiffusion(const pointVectorField& localDisp);
    void applyQualityGuardedMorph(const pointVectorField& finalPointDisp);
    pointVectorField applyRelativeGuardedMorph(const pointVectorField& finalPointDisp, bool move); // Fixed scoping
    void solveLocalLaplacian(const pointMesh& pMesh);
    scalar calculateAngleDelta(const scalar ang1, const scalar ang2);
    scalar normalizeAngle(const scalar ang);

    void wiggleSmooth(const pointMesh& pMesh);
    bool wiggleInitialized;
    //data for wiggleSmooth;
    labelHashSet wigglePoints;
    labelHashSet wigglePointsM1; //M1 inplies minus 1 aka. last step
    labelHashSet branchePoints; //All points that were soly a branch
    Map<label> pointToIdLabel;
    Map<label> pointToIdLabelM1;
    DynamicList<vector> wiggleDisplacement;
    DynamicList<scalar> wiggleLambda;
    DynamicList<scalar> wiggleLambdaM1;
    DynamicList<scalar> thetaWiggle;
    DynamicList<scalar> phiWiggle;
    DynamicList<scalar> thetaWiggleM1;
    DynamicList<scalar> phiWiggleM1;
    scalar wiggleFactor;
    scalar gradientDelta;
    scalar randomLambLimit;
    labelHashSet boundaryPoints;
    bool verboseWiggle_;
    label wiggleSurroundingLayers_;
    bool complexReward_; 
    scalar selfFactor_;

    label writeEvery_;

    scalarField calculateFaceNonOrtho() const;
    scalarField calculateFaceSkewness() const;
    scalar getMinEdgeLength(label ptI) const;

    scalar calculateCellAspectRatio(label celli) const;

    // Dynamic mesh metric functions for the Bisection Engine
    void calculateCellVolumeAndCentre
    (
        const pointField& p,
        const label celli,
        scalar& volume,
        vector& centre 
    ) const; // Fixed scoping

    scalar calculateFaceNonOrtho
    (
        const pointField& p,
        const label facei
    ) const; // Fixed scoping

    scalar calculateFaceSkewness
    (
        const pointField& p,
        const label facei
    ) const; // Fixed scoping

    // --- Additional Dynamic Mesh Metrics ---

    // Returns the face pyramid volume (positive means the face points outward correctly)
    scalar calculateFacePyramidVolume
    (
        const pointField& p,
        const label facei,
        const bool isOwner
    ) const;

    // Returns the maximum concave angle of a face in degrees (0 for convex faces)
    scalar calculateFaceConcavity
    (
        const pointField& p,
        const label facei
    ) const;

    // Returns the interpolation weight (0 to 1). Ideal is 0.5.
    scalar calculateInterpolationWeight
    (
        const pointField& p,
        const label facei
    ) const;

    // Returns the volume ratio of adjacent cells sharing the face (0 to 1). Ideal is 1.0.
    scalar calculateVolumeRatio
    (
        const pointField& p,
        const label facei
    ) const;

    // Returns the face twist as the minimum cosine between the face normal and its sub-triangles (-1 to 1). Ideal is 1.0.
    scalar calculateFaceTwist
    (
        const pointField& p,
        const label facei
    ) const;

public:

    meshFixer(fvMesh& mesh);
    void fix();
};

} // End namespace Foam

#endif