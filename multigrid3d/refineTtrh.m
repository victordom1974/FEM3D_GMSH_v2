function [Tf, info] = refineTtrh(T)
% REFINETTRH  Uniform (red) refinement of a tetrahedral mesh, Bey's algorithm.
%
%   [Tf, info] = refineTtrh(T)
%
%   T   : mesh struct as produced by importGMSH3D (fields coord, ttrh, trB,
%         domBd, domain, deg), or a FEM3Dclass object (its .mesh is used).
%         Only the vertex columns ttrh(:,1:4), trB(:,1:3) and the field deg
%         are read, so a P1 skeleton with deg = k is a valid input: the
%         output is then a proper Pk mesh (degree elevation for free).
%   Tf  : refined mesh struct with the same fields as importGMSH3D produces,
%         ready for FEM3Dclass(Tf). Every tetrahedron is split into 8 and
%         every boundary triangle into 4. Tags (domain, domBd) are inherited.
%         Pk nodes (k = T.deg, 1..4) are generated in Gmsh ordering, the
%         same used in basisNj.
%
%   info.parent      (nTtrhF x 1)  parent tetrahedron of each fine tetrahedron
%   info.child       (nTtrhF x 1)  Bey child index (1..8) of each fine tetrahedron
%   info.nodeParent  (nNodesF x 1) a coarse tetrahedron containing each fine node
%   info.nodeBaryc   (nNodesF x 4) barycentric coordinates of each fine node in
%                                  info.nodeParent, w.r.t. the vertex order of
%                                  T.ttrh(info.nodeParent, 1:4). This is all
%                                  that is needed to build the prolongation:
%                                  P(i, T.ttrh(p,j)) = Nj3D{j}(baryc(i,:)).
%
%   Notes on the algorithm (Bey 1995, Computing 55; Bey 2000, Numer. Math. 85):
%   corners are cut off at the edge midpoints and the inner octahedron is
%   split along the diagonal x13-x24. With Bey's vertex ordering of the
%   children, recursive refinement produces at most 3 similarity classes per
%   initial tetrahedron, so shape regularity is preserved on all levels.
%   Children 6 and 8 come out negatively oriented; swapping their vertices
%   1 and 3 restores detBk > 0 and keeps the diagonal rule invariant, so the
%   3-class property survives the swap (checked numerically over 5 levels).
%   Input tetrahedra must be positively oriented (detBk > 0), as Gmsh does.
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version (Bey refinement, P1-P4 nodes, boundary, info).

if isa(T, 'FEM3Dclass')
    T = T.mesh;
end
deg = T.deg;
if isempty(deg) || deg < 1 || deg > 4
    error('refineTtrh:degree', 'T.deg must be 1, 2, 3 or 4.');
end

coordC = T.coord;
nC     = size(coordC, 1);
V      = T.ttrh(:, 1:4);
nT     = size(V, 1);

% --- Orientation of the parents
dC = tetDet(V, coordC);
if any(dC <= 0)
    error('refineTtrh:orientation', ...
        '%d input tetrahedra are not positively oriented (detBk <= 0).', nnz(dC <= 0));
end

%% Edge midpoints (temporary numbering: coarse nodes first, then midpoints)
edgeLoc = [1 2; 2 3; 1 3; 1 4; 3 4; 2 4];   % Gmsh tet10 order of midpoints
E = reshape(V(:, edgeLoc.').', 2, []).';   % (6 nT) x 2, tet-major
[Eu, ~, ie] = unique(sort(E, 2), 'rows');
nE = size(Eu, 1);
mid = nC + reshape(ie, 6, nT).';           % nT x 6: m12 m23 m13 m14 m34 m24
coordTmp = [coordC; 0.5*(coordC(Eu(:,1), :) + coordC(Eu(:,2), :))];

%% Bey's eight children
x1 = V(:,1); x2 = V(:,2); x3 = V(:,3); x4 = V(:,4);
m12 = mid(:,1); m23 = mid(:,2); m13 = mid(:,3);
m14 = mid(:,4); m34 = mid(:,5); m24 = mid(:,6);

childIds = cat(3, ...
    [x1  m12 m13 m14], ...   % 1
    [m12 x2  m23 m24], ...   % 2
    [m13 m23 x3  m34], ...   % 3
    [m14 m24 m34 x4 ], ...   % 4
    [m12 m13 m14 m24], ...   % 5
    [m12 m13 m23 m24], ...   % 6
    [m13 m14 m24 m34], ...   % 7
    [m13 m23 m24 m34]);      % 8
% Barycentric coordinates (in the parent) of the children's vertices:
% same labels as above, 1..4 -> e_i, ij -> (e_i + e_j)/2
lab = {[1 12 13 14], [12 2 23 24], [13 23 3 34], [14 24 34 4], ...
       [12 13 14 24], [12 13 23 24], [13 14 24 34], [13 23 24 34]};
Bchild = zeros(4, 4, 8);
for c = 1:8
    for r = 1:4
        Bchild(r, :, c) = labelToBaryc(lab{c}(r));
    end
end
% Orientation fix: swap vertices 1 and 3 of children 6 and 8
for c = [6 8]
    childIds(:, [1 3], c) = childIds(:, [3 1], c);
    Bchild([1 3], :, c)   = Bchild([3 1], :, c);
end

Vf      = reshape(permute(childIds, [3 1 2]), 8*nT, 4);   % child-major
nTf     = 8*nT;
parent  = repelem((1:nT).', 8, 1);
child   = repmat((1:8).', nT, 1);
domainF = repelem(T.domain(:), 8, 1);

dF = tetDet(Vf, coordTmp);
if any(dF <= 0)
    error('refineTtrh:childOrientation', 'Internal error: a child is not positively oriented.');
end

%% Boundary triangles -> 4 children
hasBd = isfield(T, 'trB') && ~isempty(T.trB);
if hasBd
    Vb  = T.trB(:, 1:3);
    nB  = size(Vb, 1);
    Eb  = reshape(Vb(:, [1 2; 2 3; 1 3].').', 2, []).';
    [found, ieb] = ismember(sort(Eb, 2), Eu, 'rows');
    if ~all(found)
        error('refineTtrh:boundary', 'A boundary edge does not belong to any tetrahedron.');
    end
    midB = nC + reshape(ieb, 3, nB).';        % nB x 3: m12 m23 m13
    b1 = Vb(:,1); b2 = Vb(:,2); b3 = Vb(:,3);
    n12 = midB(:,1); n23 = midB(:,2); n13 = midB(:,3);
    bChildren = cat(3, [b1 n12 n13], [n12 b2 n23], [n13 n23 b3], [n12 n23 n13]);
    Vbf    = reshape(permute(bChildren, [3 1 2]), 4*nB, 3);
    domBdF = repelem(T.domBd(:), 4, 1);
else
    Vbf = zeros(0, 3); domBdF = zeros(0, 1);
end

%% Pk nodes on the fine mesh (Gmsh ordering), deduplicated by support + weights
L3 = localNodes3D(deg);                       % nL3 x 4 barycentric table
W3 = round(L3*deg);                           % integer weights
nL3 = size(L3, 1);

keys3 = zeros(nL3*nTf, 8);
for l = 1:nL3
    keys3(l:nL3:end, :) = nodeKeys(Vf, W3(l, :));
end
[keyU, ia, ic] = unique(keys3, 'rows');
nNodesF = size(keyU, 1);
ttrhF   = reshape(ic, nL3, nTf).';

% Coordinates: affine combination of the (temporary) vertices
idsU = keyU(:, 1:4);  wU = keyU(:, 5:8)/deg;
idsU(idsU == 0) = 1;                          % weight is 0 there anyway
coordF = wU(:,1).*coordTmp(idsU(:,1), :) + wU(:,2).*coordTmp(idsU(:,2), :) + ...
         wU(:,3).*coordTmp(idsU(:,3), :) + wU(:,4).*coordTmp(idsU(:,4), :);

% Parent tetrahedron and barycentric coordinates of each fine node
lOf = mod(ia - 1, nL3) + 1;                   % local node index
tOf = floor((ia - 1)/nL3) + 1;                % fine tetrahedron
nodeParent = parent(tOf);
nodeBaryc  = zeros(nNodesF, 4);
for c = 1:8
    sel = child(tOf) == c;
    nodeBaryc(sel, :) = L3(lOf(sel), :) * Bchild(:, :, c);
end

% Boundary Pk nodes: same keys, looked up in the volume node table
if hasBd
    L2  = localNodes2D(deg);
    W2  = round(L2*deg);
    nL2 = size(L2, 1);
    nBf = size(Vbf, 1);
    keys2 = zeros(nL2*nBf, 8);
    for l = 1:nL2
        keys2(l:nL2:end, :) = nodeKeys([Vbf, zeros(nBf, 1)], [W2(l, :), 0]);
    end
    [found, loc] = ismember(keys2, keyU, 'rows');
    if ~all(found)
        error('refineTtrh:boundaryNodes', 'Internal error: boundary node not found among volume nodes.');
    end
    trBF = reshape(loc, nL2, nBf).';
else
    trBF = zeros(0, 3);
end

%% Output mesh
Tf.nNodes = nNodesF;
Tf.deg    = deg;
Tf.coord  = coordF;
Tf.trB    = trBF;
Tf.domBd  = domBdF;
Tf.ttrh   = ttrhF;
Tf.domain = domainF;
Tf = completeMesh3D(Tf);                      % faces, normals, detBk

info.parent     = parent;
info.child      = child;
info.nodeParent = nodeParent;
info.nodeBaryc  = nodeBaryc;
info.nEdges     = nE;
end

%% ------------------------------------------------------------------------
function d = tetDet(V, coord)
v12 = coord(V(:,2), :) - coord(V(:,1), :);
v13 = coord(V(:,3), :) - coord(V(:,1), :);
v14 = coord(V(:,4), :) - coord(V(:,1), :);
d = dot(v12, cross(v13, v14, 2), 2);
end

function b = labelToBaryc(lab)
b = zeros(1, 4);
if lab < 10
    b(lab) = 1;
else
    b(floor(lab/10)) = 0.5;  b(mod(lab, 10)) = 0.5;
end
end

function K = nodeKeys(ids, w)
% One key per row of ids: vertex ids with nonzero weight, sorted, padded
% with zeros to 4 entries, followed by the weights in the same order.
n = size(ids, 1);
w = repmat(w, n, 1);
ids(w == 0) = 0;
[idsS, perm] = sort(ids, 2);
wS = w(sub2ind(size(w), repmat((1:n).', 1, 4), perm));
K = [idsS, wS];
end

function Tf = completeMesh3D(Tf)
% Derived fields, identical to the final part of importGMSH3D.
v21 = Tf.coord(Tf.trB(:,1), :) - Tf.coord(Tf.trB(:,2), :);
v31 = Tf.coord(Tf.trB(:,1), :) - Tf.coord(Tf.trB(:,3), :);
Tf.trBNormal = cross(v21, v31, 2);

nElement = size(Tf.ttrh, 1);
indAux   = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
Tf.faces = [Tf.ttrh(:, indAux(1,:)); Tf.ttrh(:, indAux(2,:)); ...
            Tf.ttrh(:, indAux(3,:)); Tf.ttrh(:, indAux(4,:))];
facesAux = sort(Tf.faces(:, 1:3), 2);
[~, p1, p3] = unique(facesAux, 'rows', 'first');
[~, p2]     = unique(facesAux, 'rows', 'last');
Tf.faces    = Tf.faces(p1, :);
ind = zeros(nElement, 4); ind(:) = 1:numel(Tf.ttrh(:, 1:4));
Tf.ttrh2faces = p3(ind);
p1 = mod(p1, nElement); p1(p1 == 0) = nElement;
p2 = mod(p2, nElement); p2(p2 == 0) = nElement;
Tf.faces2ttrh = [p1, p2];
Tf.faces2ttrh(p1 == p2, 2) = nan;

v21 = Tf.coord(Tf.faces(:,1), :) - Tf.coord(Tf.faces(:,2), :);
v31 = Tf.coord(Tf.faces(:,1), :) - Tf.coord(Tf.faces(:,3), :);
Tf.facenormal = cross(v21, v31, 2);

v12 = Tf.coord(Tf.ttrh(:,2), :) - Tf.coord(Tf.ttrh(:,1), :);
v13 = Tf.coord(Tf.ttrh(:,3), :) - Tf.coord(Tf.ttrh(:,1), :);
v14 = Tf.coord(Tf.ttrh(:,4), :) - Tf.coord(Tf.ttrh(:,1), :);
Tf.detBk = dot(v12, cross(v13, v14, 2), 2);
end

function L = localNodes3D(deg)
% Barycentric coordinates of the Pk nodes of the tetrahedron, Gmsh ordering
% (same tables as basisNj).
switch deg
    case 1
        L = eye(4);
    case 2
        L = [eye(4);
             1/2 1/2 0 0; 0 1/2 1/2 0; 1/2 0 1/2 0;
             1/2 0 0 1/2; 0 0 1/2 1/2; 0 1/2 0 1/2];
    case 3
        L = [eye(4);
             2/3 1/3 0 0; 1/3 2/3 0 0;
             0 2/3 1/3 0; 0 1/3 2/3 0;
             1/3 0 2/3 0; 2/3 0 1/3 0;
             1/3 0 0 2/3; 2/3 0 0 1/3;
             0 0 1/3 2/3; 0 0 2/3 1/3;
             0 1/3 0 2/3; 0 2/3 0 1/3;
             1/3 1/3 1/3 0; 1/3 1/3 0 1/3; 1/3 0 1/3 1/3; 0 1/3 1/3 1/3];
    case 4
        L = [eye(4);
             3/4 1/4 0 0; 2/4 2/4 0 0; 1/4 3/4 0 0;
             0 3/4 1/4 0; 0 2/4 2/4 0; 0 1/4 3/4 0;
             1/4 0 3/4 0; 2/4 0 2/4 0; 3/4 0 1/4 0;
             1/4 0 0 3/4; 2/4 0 0 2/4; 3/4 0 0 1/4;
             0 0 1/4 3/4; 0 0 2/4 2/4; 0 0 3/4 1/4;
             0 1/4 0 3/4; 0 2/4 0 2/4; 0 3/4 0 1/4;
             2/4 1/4 1/4 0; 1/4 1/4 2/4 0; 1/4 2/4 1/4 0;
             2/4 1/4 0 1/4; 1/4 2/4 0 1/4; 1/4 1/4 0 2/4;
             2/4 0 1/4 1/4; 1/4 0 1/4 2/4; 1/4 0 2/4 1/4;
             0 1/4 1/4 2/4; 0 2/4 1/4 1/4; 0 1/4 2/4 1/4;
             1/4 1/4 1/4 1/4];
end
end

function L = localNodes2D(deg)
% Barycentric coordinates of the Pk nodes of the triangle, Gmsh ordering
% (consistent with the Ns2D functions in basisNj).
switch deg
    case 1
        L = eye(3);
    case 2
        L = [eye(3); 1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2];
    case 3
        L = [eye(3);
             2/3 1/3 0; 1/3 2/3 0; 0 2/3 1/3; 0 1/3 2/3; 1/3 0 2/3; 2/3 0 1/3;
             1/3 1/3 1/3];
    case 4
        L = [eye(3);
             3/4 1/4 0; 2/4 2/4 0; 1/4 3/4 0;
             0 3/4 1/4; 0 2/4 2/4; 0 1/4 3/4;
             1/4 0 3/4; 2/4 0 2/4; 3/4 0 1/4;
             2/4 1/4 1/4; 1/4 2/4 1/4; 1/4 1/4 2/4];
end
end
