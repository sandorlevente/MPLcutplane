foamToVTK -faceSet nonOrthoFaces
foamToVTK -faceSet highAspectRatioCells
foamToVTK -faceSet zeroVolumeCells
foamToVTK -faceSet wrongOrientedFaces


surfaceTransformPoints -rollPitchYaw "(90 0 0)"
surfaceTransformPoints -rotate-angle "((1 0 0) 90)"

transformPoints -scale '(1.33333 1.33333 1.3333)'
transformPoints -rollPitchYaw "(90 0 0)"