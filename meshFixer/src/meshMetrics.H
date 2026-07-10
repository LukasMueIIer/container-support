// ------------------------------------------------------------------------- //
// Helper: Dynamically calculate a cell's center and volume
// ------------------------------------------------------------------------- //
void meshFixer::calculateCellVolumeAndCentre
(
    const pointField& p,
    const label celli,
    scalar& volume,
    vector& centre
) const
{
    const cell& cFaces = mesh_.cells()[celli];
    
    // 1. Estimate center via face centers
    vector cEst = Zero;
    for (const label facei : cFaces)
    {
        cEst += mesh_.faces()[facei].centre(p);
    }
    cEst /= cFaces.size();

    // 2. Sum face-pyramid contributions
    centre = Zero;
    volume = 0.0;

    for (const label facei : cFaces)
    {
        const face& f = mesh_.faces()[facei];
        vector fc = f.centre(p);
        vector fa = f.areaNormal(p);
        
        scalar pyr3Vol = fa & (fc - cEst);
        if (mesh_.faceOwner()[facei] != celli)
        {
            pyr3Vol = -pyr3Vol;
        }

        volume += pyr3Vol;
        centre += pyr3Vol * (0.75 * fc + 0.25 * cEst);
    }

    if (mag(volume) > VSMALL)
    {
        centre /= volume;
    }
    else
    {
        centre = cEst;
    }
    volume /= 3.0;
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Non-Orthogonality
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceNonOrtho
(
    const pointField& p,
    const label facei
) const
{
    // 1. OPENFOAM COMPLIANCE: Solid walls have perfect orthogonality (0 degrees)
    if (!mesh_.isInternalFace(facei))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.whichPatch(facei);
        
        // If it is an uncoupled boundary (wall, patch, empty), OpenFOAM ignores it.
        if (!pbm[patchi].coupled())
        {
            return 0.0; 
        }
    }

    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    
    if (mag(s) < VSMALL) return GREAT; 

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);

    vector d;
    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        // For coupled patches (MPI/Cyclic), we approximate using the face center.
        // To be 100% compliant, you would need to sync dynamic cell centers 
        // across processors here, but for a local smoother, this approximation 
        // is standard to avoid massive MPI overhead inside the loop.
        d = f.centre(p) - ownCc;
    }

    // 2. OPENFOAM COMPLIANCE: Calculate Cosine and Clamp
    scalar dDotS = (d & s) / (mag(d) * mag(s) + VSMALL);
    dDotS = max(-1.0, min(1.0, dDotS)); 
    
    // Return the Angle (to match your qBad_ logic and checkMesh output)
    return radToDeg(acos(dDotS));
}

// ------------------------------------------------------------------------- //
// Helper: Calculate Skewness
// ------------------------------------------------------------------------- //
scalar meshFixer::calculateFaceSkewness
(
    const pointField& p,
    const label facei
) const
{
    const face& f = mesh_.faces()[facei];
    vector s = f.areaNormal(p);
    if (mag(s) < VSMALL) return GREAT;

    const label own = mesh_.faceOwner()[facei];
    scalar ownVol; vector ownCc;
    calculateCellVolumeAndCentre(p, own, ownVol, ownCc);
    
    vector fc = f.centre(p);
    vector d;
    vector dOwn = fc - ownCc;

    if (mesh_.isInternalFace(facei))
    {
        const label nei = mesh_.faceNeighbour()[facei];
        scalar neiVol; vector neiCc;
        calculateCellVolumeAndCentre(p, nei, neiVol, neiCc);
        d = neiCc - ownCc;
    }
    else
    {
        d = fc - ownCc;
    }

    // Intersect d with the face plane
    scalar t = (s & dOwn) / ((s & d) + VSMALL);
    vector pInt = ownCc + t * d;

    return mag(fc - pInt) / (mag(d) + VSMALL);
}