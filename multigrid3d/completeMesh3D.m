function T = completeMesh3D(T)
% COMPLETEMESH3D  Fill in the derived fields of a mesh struct.
%
%   T = completeMesh3D(T) takes a mesh struct with at least coord, ttrh,
%   trB, domBd, domain and deg, and computes nNodes, trBNormal, faces,
%   ttrh2faces, faces2ttrh, facenormal and detBk exactly as importGMSH3D
%   does, so that FEM3Dclass(T) works. Nodes not referenced by ttrh
%   (ghost nodes, e.g. the P2 edge nodes after taking the P1 skeleton) are
%   removed and the remaining ones renumbered. If the field refinement is
%   missing, an empty level-0 record is added (coarsest mesh of a hierarchy).
%
%   The code of the derived fields is the final part of importGMSH3D,
%   reproduced here so that meshes not read from Gmsh (refined, hand-made)
%   get the same structure.
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version (extracted from refineTtrh / importGMSH3D).
%   2026-10-08  Ghost nodes removed and renumbered.

if ~isfield(T, 'trB') || isempty(T.trB)
    T.trB = zeros(0, 3);  T.domBd = zeros(0, 1);
end
% Ghost nodes: keep only nodes referenced by ttrh, renumber
used = unique(T.ttrh(:));
if numel(used) ~= size(T.coord, 1) || any(used(:).' ~= 1:numel(used))
    map = zeros(max(size(T.coord, 1), max(used)), 1);
    map(used) = 1:numel(used);
    if ~isempty(T.trB) && any(map(T.trB(:)) == 0)
        error('completeMesh3D:boundary', 'trB references nodes not in ttrh.');
    end
    T.ttrh  = reshape(map(T.ttrh), size(T.ttrh));
    T.trB   = reshape(map(T.trB), size(T.trB));
    T.coord = T.coord(used, :);
end
T.nNodes = numel(used);

v21 = T.coord(T.trB(:,1), :) - T.coord(T.trB(:,2), :);
v31 = T.coord(T.trB(:,1), :) - T.coord(T.trB(:,3), :);
T.trBNormal = cross(v21, v31, 2);

nElement = size(T.ttrh, 1);
indAux   = [2 3 4; 1 4 3; 1 2 4; 1 3 2];
T.faces  = [T.ttrh(:, indAux(1,:)); T.ttrh(:, indAux(2,:)); ...
            T.ttrh(:, indAux(3,:)); T.ttrh(:, indAux(4,:))];
facesAux = sort(T.faces(:, 1:3), 2);
[~, p1, p3] = unique(facesAux, 'rows', 'first');
[~, p2]     = unique(facesAux, 'rows', 'last');
T.faces     = T.faces(p1, :);
ind = zeros(nElement, 4); ind(:) = 1:numel(T.ttrh(:, 1:4));
T.ttrh2faces = p3(ind);
p1 = mod(p1, nElement); p1(p1 == 0) = nElement;
p2 = mod(p2, nElement); p2(p2 == 0) = nElement;
T.faces2ttrh = [p1, p2];
T.faces2ttrh(p1 == p2, 2) = nan;

v21 = T.coord(T.faces(:,1), :) - T.coord(T.faces(:,2), :);
v31 = T.coord(T.faces(:,1), :) - T.coord(T.faces(:,3), :);
T.facenormal = cross(v21, v31, 2);

v12 = T.coord(T.ttrh(:,2), :) - T.coord(T.ttrh(:,1), :);
v13 = T.coord(T.ttrh(:,3), :) - T.coord(T.ttrh(:,1), :);
v14 = T.coord(T.ttrh(:,4), :) - T.coord(T.ttrh(:,1), :);
T.detBk = dot(v12, cross(v13, v14, 2), 2);

if ~isfield(T, 'refinement')
    T.refinement = struct('level', 0, 'parent', [], 'child', [], ...
                          'nodeParent', [], 'nodeBaryc', [], 'nEdges', []);
end
end
