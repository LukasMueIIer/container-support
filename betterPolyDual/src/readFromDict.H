IOdictionary advPolyDualDict
(
    IOobject
    (
        "advPolyDualDict", // Name of the dictionary file
        mesh.time().system(),   // Location (looks in the 'system' folder)
        mesh,                   // Object registry
        IOobject::MUST_READ,    // Throw an error if the file is missing
        IOobject::NO_WRITE      // We only want to read, not modify the file
    )
);

Info << "Reading adv(anced)PolyDualDict\n" << endl;

// 1. Declare your OpenFOAM integer (label)
label flagBoundaryHandlingMethod = 0; // Defaulting to 0 (none)

// 2. Read the boundary handling method as a word
// Using "none" as a safe fallback if the keyword is completely missing
word bhm = advPolyDualDict.lookupOrDefault<word>("boundaryHandlingMethod", "none");

// 3. Map the word to the corresponding integer flag
if (bhm == "none")
{
    flagBoundaryHandlingMethod = 0;
    Info << "Boundary handling method: none (0)" << endl;
}
else if (bhm == "keepAll")
{
    flagBoundaryHandlingMethod = 1;
    Info << "Boundary handling method: keepAll (1)" << endl;
}
else if (bhm == "cellQualityBased")
{
    flagBoundaryHandlingMethod = 2;
    Info << "Boundary handling method: cellQualityBased (2)" << endl;
}
else
{
    // 4. Robust error handling for invalid user input
    WarningInFunction 
        << "Unknown boundaryHandelingMethod: '" << bhm 
        << "'. Valid options are: none, keepAll, cellQualityBased." << nl
        << "Defaulting to 'none' (0)." << endl;
        
    flagBoundaryHandlingMethod = 0;
}

// 23. Read the flag for extruded climb boundary splitting
// Cast to <Switch> to allow true/false, yes/no, on/off, then store as a standard C++ bool
bool extrudedClimbBoundarySplitting = advPolyDualDict.lookupOrDefault<Switch>("extrudedClimbBoundarySplitting", false);

// Output to the console to confirm
Info << "Extruded climb boundary splitting: " << extrudedClimbBoundarySplitting << endl;

// ------------------------------------------------------------------------- //

// 24. Read the maximum number of climbing iterations
// Read directly as an OpenFOAM label (integer), with 10 as your fallback safeguard
label maxClimbingIterations = advPolyDualDict.lookupOrDefault<label>("maxClimbingIterations", 10);

// Output to the console to confirm
Info << "Max climbing iterations: " << maxClimbingIterations << endl;

scalar boundaryOrthoDetectionLimit = advPolyDualDict.lookupOrDefault<scalar>("boundaryOrthoDetectionLimit", 20.0);

bool detectionVerbose = advPolyDualDict.lookupOrDefault<Switch>("detectionVerbose", false);

// Output to the console to confirm it was read correctly
// Note: 'detectionVerbose' will print as 1 (true) or 0 (false) in the terminal
Info << "Verbose detection output: " << detectionVerbose << endl;


// 7. Declare the integer flag for non-tet cell handling
label flagNonTetCellMethod = 0; // Defaulting to 0 (none)

// 8. Read the method as a word from the dictionary
// Using "none" as a safe fallback
word ntcm = advPolyDualDict.lookupOrDefault<word>("nonTetCellMethod", "none");

// 9. Map the word to the corresponding integer flag
if (ntcm == "none")
{
    flagNonTetCellMethod = 0;
    Info << "Non-tet cell handling method: none (0)" << endl;
}
else if (ntcm == "singleSplit")
{
    flagNonTetCellMethod = 1;
    Info << "Non-tet cell handling method: singleSplit (1)" << endl;
}
else if (ntcm == "multiSplit")
{
    flagNonTetCellMethod = 2;
    Info << "Non-tet cell handling method: multiSplit (2)" << endl;
}
else
{
    // 10. Robust error handling for invalid user input
    WarningInFunction 
        << "Unknown nonTetCellMethod: '" << ntcm 
        << "'. Valid options are: none, singleSplit, multiSplit." << nl
        << "Defaulting to 'none' (0)." << endl;
        
    flagNonTetCellMethod = 0;
}


// 11. Declare the integer flag for interface treatment method
label flagInterfaceTreatmentMethod = 0; // Defaulting to 0 (none)

// 12. Read the interface treatment method as a word
word itm = advPolyDualDict.lookupOrDefault<word>("interfaceTreatmentMethod", "none");

// 13. Map the word to the corresponding integer flag
if (itm == "none")
{
    flagInterfaceTreatmentMethod = 0;
    Info << "Interface treatment method: none (0)" << endl;
}
else if (itm == "singleSplit")
{
    flagInterfaceTreatmentMethod = 1;
    Info << "Interface treatment method: singleSplit (1)" << endl;
}
else if (itm == "multiSplit")
{
    flagInterfaceTreatmentMethod = 2;
    Info << "Interface treatment method: multiSplit (2)" << endl;
}
else
{
    // Robust error handling for invalid user input
    WarningInFunction 
        << "Unknown interfaceTreatmentMethod: '" << itm 
        << "'. Valid options are: none, singleSplit, multiSplit." << nl
        << "Defaulting to 'none' (0)." << endl;
        
    flagInterfaceTreatmentMethod = 0;
}

// ------------------------------------------------------------------------- //

// 14. Read the interface buffer level
// We read this directly as an OpenFOAM label (integer). 
// Defaulting to 0 disables it, matching your documentation.
label interfaceBufferLevel = advPolyDualDict.lookupOrDefault<label>("interfaceBufferLevel", 0);

// Output to the console to confirm it was read correctly
Info << "Interface buffer level: " << interfaceBufferLevel << endl;


// 15. Declare the integer flag for hex cell treatment
label flagHexCellMethod = 0; // Defaulting to 0 (none)

// 16. Read the hex cell method as a word
word hcm = advPolyDualDict.lookupOrDefault<word>("hexCellMethod", "none");

// 17. Map the word to the corresponding integer flag
if (hcm == "none")
{
    flagHexCellMethod = 0;
    Info << "Hex cell handling method: none (0)" << endl;
}
else if (hcm == "singleSplit")
{
    flagHexCellMethod = 1;
    Info << "Hex cell handling method: singleSplit (1)" << endl;
}
else if (hcm == "multiSplit")
{
    flagHexCellMethod = 2;
    Info << "Hex cell handling method: multiSplit (2)" << endl;
}
else
{
    // 18. Robust error handling for invalid user input
    WarningInFunction 
        << "Unknown hexCellMethod: '" << hcm 
        << "'. Valid options are: none, singleSplit, multiSplit." << nl
        << "Defaulting to 'none' (0)." << endl;
        
    flagHexCellMethod = 0;
}

// 19. Declare the integer flag for hex interface treatment
label flagHexInterfaceTreatmentMethod = 0; // Defaulting to 0 (none)

// 20. Read the hex interface treatment method as a word
word hitm = advPolyDualDict.lookupOrDefault<word>("hexInterfaceTreatmentMethod", "none");

// 21. Map the word to the corresponding integer flag
if (hitm == "none")
{
    flagHexInterfaceTreatmentMethod = 0;
    Info << "Hex interface treatment method: none (0)" << endl;
}
else if (hitm == "singleSplit")
{
    flagHexInterfaceTreatmentMethod = 1;
    Info << "Hex interface treatment method: singleSplit (1)" << endl;
}
else if (hitm == "multiSplit")
{
    flagHexInterfaceTreatmentMethod = 2;
    Info << "Hex interface treatment method: multiSplit (2)" << endl;
}
else
{
    // 22. Robust error handling for invalid user input
    WarningInFunction 
        << "Unknown hexInterfaceTreatmentMethod: '" << hitm 
        << "'. Valid options are: none, singleSplit, multiSplit." << nl
        << "Defaulting to 'none' (0)." << endl;
        
    flagHexInterfaceTreatmentMethod = 0;
}

// 25. Read the flag for quality-based cell splitting
bool useQualitySplitting = advPolyDualDict.lookupOrDefault<Switch>("useQualitySplitting", true);

Info << "Use quality splitting: " << useQualitySplitting << endl;

// 26. Read the quality limits
// These are numerical thresholds, so we read them as 'scalar' (double precision)
scalar skewnessLimit      = advPolyDualDict.lookupOrDefault<scalar>("skewnessLimit", 0.4);
scalar nonOrthoAngleLimit = advPolyDualDict.lookupOrDefault<scalar>("nonOrthoAngleLimit", 45.0);
scalar minTetVolumeLimit  = advPolyDualDict.lookupOrDefault<scalar>("minTetVolumeLimit", 1e-4);
scalar aspectRatioLimit   = advPolyDualDict.lookupOrDefault<scalar>("aspectRatioLimit", 20.0);

// Output to the console only if quality splitting is turned on
if (useQualitySplitting)
{
    Info << "  - Skewness limit        : " << skewnessLimit << nl
         << "  - Non-ortho angle limit : " << nonOrthoAngleLimit << nl
         << "  - Min tet volume limit  : " << minTetVolumeLimit << nl
         << "  - Aspect ratio limit    : " << aspectRatioLimit << endl;
}

// 27. Read the quality transition levels (number of buffer cell layers)
// We read this as a label (integer). Defaulting to 1 based on your notes.
label qualityTransitionLevels = advPolyDualDict.lookupOrDefault<label>("qualityTransitionLevels", 1);

// Optional: Output to the console. 
// It is good practice to put this inside your previous 'if (useQualitySplitting)' 
// block so it only prints when quality splitting is actually active!
if (useQualitySplitting)
{
    Info << "  - Quality transition levels : " << qualityTransitionLevels << endl;
}
