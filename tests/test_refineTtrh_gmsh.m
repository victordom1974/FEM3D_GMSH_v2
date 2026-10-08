% TEST_REFINETTRH_GMSH
% Refines a real Gmsh mesh: Examples/Elasticity/beam.msh, a P2 mesh of the
% box [10,30] x [1.5,4.5] x [0,150] with physical surface tags.
% Checks volume (9000), boundary area (7020), inheritance of the boundary
% tags (each tag count x4 per level), node count, mass matrix total,
% and exact reproduction of quadratics by the nested prolongation.
% Run from the repository root with addpath(genpath(pwd)).
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version.

clear; clc;

meshFile = fullfile('Examples', 'Elasticity', 'beam.msh');
obj0 = FEM3Dclass(meshFile);
T0 = obj0.mesh;
fprintf('beam.msh: deg %d, %d ttrh, %d nodes, %d boundary triangles, tags %s\n', ...
    T0.deg, size(T0.ttrh,1), T0.nNodes, size(T0.trB,1), mat2str(unique(T0.domBd)'));

volExact  = 20*3*150;
areaExact = 2*(20*3 + 20*150 + 3*150);
report(T0, 0, volExact, areaExact);

T = T0;
for lev = 1:2
    tic;
    [T, info] = refineTtrh(T);
    t = toc;
    report(T, lev, volExact, areaExact);
    fprintf('          refinement time %.2f s\n', t);
    % boundary tags inherited
    tags = unique(T0.domBd);
    for k = 1:numel(tags)
        assert(nnz(T.domBd == tags(k)) == 4^lev*nnz(T0.domBd == tags(k)), 'domBd inheritance');
    end
end

% Assembly on each level: total mass must be the volume
obj = FEM3Dclass(T);
M = obj.femMassMatrix(@(x,y,z) 1 + 0*x);
fprintf('finest level: sum(M) - volume = %.2e\n', sum(M(:)) - volExact);

% Prolongation level 1 -> level 2 reproduces quadratics (deg = 2)
[T1, ~]     = refineTtrh(T0);
[T2, info2] = refineTtrh(T1);
obj1 = FEM3Dclass(T1);
P = prolongationFromInfo(obj1, T2, info2);
f = @(x) 1 + 2*x(:,1) - 3*x(:,2) + 0.5*x(:,3) + x(:,1).^2 - x(:,2).*x(:,3) + x(:,3).^2;
err = norm(P*f(T1.coord) - f(T2.coord), inf) / norm(f(T2.coord), inf);
fprintf('prolongation %d x %d: relative reproduction error %.2e\n', size(P,1), size(P,2), err);

% ------------------------------------------------------------------------
function report(T, lev, volExact, areaExact)
V = T.ttrh(:,1:4);
nV = numel(unique(V(:)));
E = reshape(V(:, [1 2; 2 3; 1 3; 1 4; 3 4; 2 4].').', 2, []).';
nE = size(unique(sort(E,2), 'rows'), 1);
fprintf('level %d: %7d ttrh, %8d nodes (nV+nE = %8d), min detBk %.2e, vol-err %.2e, area-err %.2e\n', ...
    lev, size(V,1), T.nNodes, nV + nE, min(T.detBk), ...
    sum(T.detBk)/6 - volExact, sum(vecnorm(T.trBNormal,2,2))/2 - areaExact);
end

function P = prolongationFromInfo(objC, Tf, info)
Tc = objC.mesh; nL = size(Tc.ttrh,2); nF = Tf.nNodes; b = info.nodeBaryc;
vals = zeros(nF, nL);
for j = 1:nL
    vals(:,j) = objC.Nj3D{j}(b(:,1), b(:,2), b(:,3), b(:,4));
end
cols = Tc.ttrh(info.nodeParent, :);
rows = repmat((1:nF).', 1, nL);
keep = abs(vals) > 1e-14;
P = sparse(rows(keep), cols(keep), vals(keep), nF, Tc.nNodes);
end
