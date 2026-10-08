# Branch `multigrid3d`

Geometric multigrid for the 3D finite element package in `@FEM3Dclass`
(tetrahedral Lagrange elements P1-P4, meshes from Gmsh 4.1).

Fork of [aleduques/FEM3D_GMSH](https://github.com/aleduques/FEM3D_GMSH),
the original package by Alejandro Duque Salazar (bachelor's thesis, UPNA).
The class folder `@FEM3Dclass/` is kept as in the original, apart from the
minimal changes listed below. All new code lives outside it.

## Layout

    @FEM3Dclass/      original package
    multigrid3d/      refineTtrh.m         uniform (red) refinement, Bey's algorithm, P1-P4
                      completeMesh3D.m     derived mesh fields (faces, normals, detBk)
                      prolongation3D.m     nested prolongation from refineTtrh info   (to do)
                      @MG3Dsolver/         multigrid V-cycle / preconditioner            (to do)
    tests/            test_refineTtrh.m        Kuhn cube, no Gmsh needed
                      test_refineTtrh_gmsh.m   Examples/Elasticity/beam.msh (P2 box)
    Examples/         original examples

## Changes to `@FEM3Dclass` (planned, minimal)

- `basisNj`: also return the table of barycentric node coordinates.
- `importGMSH3D`: move the computation of derived fields (faces, normals,
  detBk) to a separate function so that refined meshes can reuse it.

## Usage

    addpath(genpath(pwd))        % from the repository root
    run tests/test_refineTtrh.m

Code and comments in English. Each file carries its last-modification date
and a changelog in the header.
