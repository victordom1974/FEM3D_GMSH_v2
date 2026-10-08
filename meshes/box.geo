// box.geo - parameterized box for the multigrid tests.
//
// Usage (from gmshBox.m or the command line):
//   gmsh -3 -order 2 -format msh41 -setnumber Lx 20 -setnumber Ly 3 -setnumber Lz 150 \
//        -setnumber lc 3 box.geo -o box.msh
//
// Physical tags: surfaces xmin, xmax, ymin, ymax, zmin, zmax; volume "box".
// Mesh optimization is switched on so that no slivers survive: the quality of
// the coarse mesh is inherited by every level of the refinement hierarchy.
//
// Last modified: 2026-10-08
// Changelog:
//   2026-10-08  First version.

SetFactory("OpenCASCADE");

DefineConstant[
  x0 = {0,   Name "Origin x"},
  y0 = {0,   Name "Origin y"},
  z0 = {0,   Name "Origin z"},
  Lx = {1,   Name "Length x"},
  Ly = {1,   Name "Length y"},
  Lz = {1,   Name "Length z"},
  lc = {0.25, Name "Element size"}
];

Box(1) = {x0, y0, z0, Lx, Ly, Lz};

eps = 1e-6 * (Lx + Ly + Lz);
sxmin() = Surface In BoundingBox{x0-eps, y0-eps, z0-eps, x0+eps,    y0+Ly+eps, z0+Lz+eps};
sxmax() = Surface In BoundingBox{x0+Lx-eps, y0-eps, z0-eps, x0+Lx+eps, y0+Ly+eps, z0+Lz+eps};
symin() = Surface In BoundingBox{x0-eps, y0-eps, z0-eps, x0+Lx+eps, y0+eps,    z0+Lz+eps};
symax() = Surface In BoundingBox{x0-eps, y0+Ly-eps, z0-eps, x0+Lx+eps, y0+Ly+eps, z0+Lz+eps};
szmin() = Surface In BoundingBox{x0-eps, y0-eps, z0-eps, x0+Lx+eps, y0+Ly+eps, z0+eps};
szmax() = Surface In BoundingBox{x0-eps, y0-eps, z0+Lz-eps, x0+Lx+eps, y0+Ly+eps, z0+Lz+eps};

Physical Surface("xmin", 1) = sxmin();
Physical Surface("xmax", 2) = sxmax();
Physical Surface("ymin", 3) = symin();
Physical Surface("ymax", 4) = symax();
Physical Surface("zmin", 5) = szmin();
Physical Surface("zmax", 6) = szmax();
Physical Volume("box", 1) = {1};

Mesh.MeshSizeMin = lc;
Mesh.MeshSizeMax = lc;
Mesh.Algorithm3D = 1;        // Delaunay
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.OptimizeThreshold = 0.3;
Mesh.MshFileVersion = 4.1;
