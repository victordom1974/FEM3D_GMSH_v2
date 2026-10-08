% TEST_REFINETTRH
% Starts from the unit cube split into 6 Kuhn tetrahedra (no Gmsh needed),
% refines it several times for deg = 1..4 and checks:
%   - all detBk > 0, total volume 1, boundary area 6
%   - node count (deg 1: nV; deg 2: nV + nE)
%   - shape regularity: min quality constant across levels (Bey's 3 classes)
%   - FEM3Dclass assembly: sum(M) = volume, S*ones = 0
%   - nested prolongation built from info reproduces P_deg polynomials exactly
% Requires FEM3Dclass and its methods on the path.
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version.
%   2026-10-08  Use completeMesh3D and T.refinement.
%   2026-10-08  Use prolongation3D.
%   2026-10-08  Use meshQuality3D.

clear; clc;

for deg = 1:4
    fprintf('\n=== deg = %d ===\n', deg);
    T0 = cubeKuhnMesh(deg);
    T  = T0;
    q0 = meshQuality3D(T);
    nLev = 3;
    for lev = 1:nLev
        T = refineTtrh(T);
        q = meshQuality3D(T);
        area = sum(vecnorm(T.trBNormal, 2, 2))/2;
        fprintf('level %d: %7d ttrh, %8d nodes, min detBk %.2e, vol %.12f, area %.12f, min quality %.4f (coarse %.4f)\n', ...
            lev, size(T.ttrh,1), T.nNodes, min(T.detBk), sum(T.detBk)/6, area, min(q), min(q0));
        % Node count: nV + (deg-1) nE + (deg-1)(deg-2)/2 nF + (deg-1)(deg-2)(deg-3)/6 nT
        V = T.ttrh(:, 1:4);
        nV = numel(unique(V(:)));
        E = reshape(V(:, [1 2; 2 3; 1 3; 1 4; 3 4; 2 4].').', 2, []).';
        nE = size(unique(sort(E, 2), 'rows'), 1);
        nFc = size(T.faces, 1);
        nTt = size(V, 1);
        expected = nV + (deg-1)*nE + (deg-1)*(deg-2)/2*nFc + (deg-1)*(deg-2)*(deg-3)/6*nTt;
        assert(T.nNodes == expected && max(T.ttrh(:)) == T.nNodes, 'node count mismatch');
    end

    % Assembly on the finest level
    obj = FEM3Dclass(T);
    M = obj.femMassMatrix(@(x,y,z) 1 + 0*x);
    S = obj.femStressMatrix();
    e = ones(T.nNodes, 1);
    fprintf('assembly: sum(M) - 1 = %.2e, |S*1| = %.2e, symmetry |S-S''| = %.2e\n', ...
        sum(M(:)) - 1, norm(S*e, inf), norm(S - S', inf));

    % Prolongation coarse (two levels up) -> fine, polynomial reproduction
    Tc = refineTtrh(T0);
    Tf = refineTtrh(Tc);
    objC = FEM3Dclass(Tc);
    P = prolongation3D(objC, Tf);
    f = @(x) (1 + 2*x(:,1) - 3*x(:,2) + 0.5*x(:,3)).^deg;     % exact in P_deg
    err = norm(P*f(Tc.coord) - f(Tf.coord), inf);
    fprintf('prolongation: size %d x %d, nnz/row %.2f, polynomial reproduction error %.2e\n', ...
        size(P,1), size(P,2), nnz(P)/size(P,1), err);
end

% ------------------------------------------------------------------------
cols = Tc.ttrh(info.nodeParent, :);
rows = repmat((1:nF).', 1, nL);
keep = abs(vals) > 1e-14;
P = sparse(rows(keep), cols(keep), vals(keep), nF, Tc.nNodes);
end

q = 6*sqrt(2)*(T.detBk/6) ./ max(h, [], 2).^3;
end

function T = cubeKuhnMesh(deg)
% Unit cube, 6 Kuhn tetrahedra, P1 skeleton with T.deg = deg (refineTtrh
% generates the Pk nodes). Derived fields via completeMesh3D.
[x, y, z] = ndgrid([0 1]);
coord = [x(:) y(:) z(:)];
id = @(p) 1 + p(1) + 2*p(2) + 4*p(3);
perms3 = perms(1:3);
ttrh = zeros(6, 4);
for k = 1:6
    p = [0 0 0]; v = id(p);
    for s = 1:3
        p(perms3(k, s)) = 1; v(end+1) = id(p); %#ok<AGROW>
    end
    ttrh(k, :) = v;
end
d = tetDetLocal(ttrh, coord);
ttrh(d < 0, [1 3]) = ttrh(d < 0, [3 1]);
% Boundary triangles: faces belonging to a single tetrahedron
indAux = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
faces = [ttrh(:, indAux(1,:)); ttrh(:, indAux(2,:)); ttrh(:, indAux(3,:)); ttrh(:, indAux(4,:))];
[~, ~, j] = unique(sort(faces, 2), 'rows');
cnt = accumarray(j, 1);
trB = faces(cnt(j) == 1, :);

T.deg = deg; T.coord = coord;
T.ttrh = ttrh; T.domain = ones(6, 1);
T.trB = trB;   T.domBd = ones(size(trB, 1), 1);
T = completeMesh3D(T);
end

function d = tetDetLocal(V, coord)
v12 = coord(V(:,2), :) - coord(V(:,1), :);
v13 = coord(V(:,3), :) - coord(V(:,1), :);
v14 = coord(V(:,4), :) - coord(V(:,1), :);
d = dot(v12, cross(v13, v14, 2), 2);
end
