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
                      prolongation3D.m     nested prolongation between consecutive levels
                      MG3D.m               multigrid class: V-cycle, smoothers (sgs/gs/jacobi/ilu), pcg/gmres
                      gmshBox.m            box mesh via Gmsh (meshes/box.geo), optimized
                      meshQuality3D.m      per-tetrahedron shape quality
    meshes/           box.geo              parameterized box with physical tags
    tests/            test_refineTtrh.m        Kuhn cube, no Gmsh needed
                      test_refineTtrh_gmsh.m   Examples/Elasticity/beam.msh (P2 box)
                      test_MG3D.m              -Lap u + u = f, Neumann, Gmsh box (fallback beam.msh)
    Examples/         original examples

## Changes to `@FEM3Dclass` (planned, minimal)

- `basisNj`: also return the table of barycentric node coordinates.
- `importGMSH3D`: move the computation of derived fields (faces, normals,
  detBk) to a separate function so that refined meshes can reuse it.

## Usage

    addpath(genpath(pwd))        % from the repository root
    run tests/test_refineTtrh.m

Mesh quality matters: Bey's refinement preserves the similarity classes of
the coarse mesh, so slivers in the coarse mesh are reproduced on every
level and ruin point smoothers (beam.msh: rho ~ 0.8). Use an optimized
coarse mesh (gmshBox) and check meshQuality3D.

Code and comments in English. Each file carries its last-modification date
and a changelog in the header.
