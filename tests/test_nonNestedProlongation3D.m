% TEST_NONNESTEDPROLONGATION3D
% Builds two unrelated tet meshes of the unit cube (no Gmsh needed), lifts
% them to P2 (Gmsh tet10 ordering) and checks that the non-nested prolongation
% reproduces polynomials exactly:
%   P1 -> P1 : linear functions
%   P2 -> P2 : quadratic functions
%   P1 -> P2 : linear functions (p-multigrid transfer)
% Expected errors: O(1e-13).
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version.

clear; clc;
rng(1);

coarse1 = boxTetMesh(4);   % 4 points per side
fine1   = boxTetMesh(9);   % 9 points per side, unrelated to the coarse one

coarse2 = p1top2(coarse1);
fine2   = p1top2(fine1);

flin  = @(x) 1 + 2*x(:,1) - 3*x(:,2) + 0.5*x(:,3);
fquad = @(x) flin(x) + x(:,1).^2 + 2*x(:,1).*x(:,2) - x(:,3).^2 + x(:,2).*x(:,3);

checkTransfer('P1 -> P1, linear',    coarse1, fine1, flin);
checkTransfer('P2 -> P2, quadratic', coarse2, fine2, fquad);
checkTransfer('P1 -> P2, linear',    coarse1, fine2, flin);

% ------------------------------------------------------------------------
function checkTransfer(label, coarse, fine, f)
[P, info] = nonNestedProlongation3D(coarse, fine);
usedF = unique(fine.elems(:));
uc = f(coarse.nodes);
uf = f(fine.nodes);
err = max(abs(P*uc - uf) ./ (1 + abs(uf)));
err = max(err(usedF));
fprintf('%-24s size %6d x %6d, nnz/row %5.2f, outside %d, viol %.1e, err %.2e\n', ...
    label, size(P,1), size(P,2), nnz(P)/numel(usedF), info.nOutside, info.maxViolation, err);
end

function m = boxTetMesh(n)
% Unit cube, n^3 grid points, interior ones jittered to avoid degenerate
% Delaunay configurations; boundary points stay exactly on the faces.
[x, y, z] = ndgrid(linspace(0, 1, n));
pts = [x(:), y(:), z(:)];
h = 1/(n-1);
inner = all(pts > 0 & pts < 1, 2);
pts(inner, :) = pts(inner, :) + 0.2*h*(rand(nnz(inner), 3) - 0.5);
m.nodes = pts;
m.elems = delaunayn(pts, {'Qt', 'Qbb', 'Qc', 'Qz'});
end

function m2 = p1top2(m1)
% Adds edge midpoints and builds Gmsh tet10 connectivity.
T = m1.elems;
nV = size(m1.nodes, 1);
edges = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
E = zeros(size(T,1)*6, 2);
for k = 1:6
    E((k-1)*size(T,1)+(1:size(T,1)), :) = T(:, edges(k,:));
end
[ue, ~, j] = unique(sort(E, 2), 'rows');
mid = 0.5*(m1.nodes(ue(:,1), :) + m1.nodes(ue(:,2), :));
m2.nodes = [m1.nodes; mid];
m2.elems = [T, nV + reshape(j, size(T,1), 6)];
end
