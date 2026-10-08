function q = meshQuality3D(T)
% MESHQUALITY3D  Shape quality of each tetrahedron.
%
%   q = meshQuality3D(T) returns 6*sqrt(2)*volume / (longest edge)^3 for
%   each tetrahedron of the mesh struct (or FEM3Dclass) T: 1 for the
%   regular tetrahedron, 0 for a degenerate one. Values below ~0.1 are
%   slivers, which point smoothers cannot handle.
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version (moved from the tests).

if isa(T, 'FEM3Dclass'), T = T.mesh; end
c = T.coord; V = T.ttrh(:, 1:4);
v12 = c(V(:,2), :) - c(V(:,1), :);
v13 = c(V(:,3), :) - c(V(:,1), :);
v14 = c(V(:,4), :) - c(V(:,1), :);
vol = abs(dot(v12, cross(v13, v14, 2), 2)) / 6;
e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
h = zeros(size(V, 1), 6);
for k = 1:6
    h(:, k) = vecnorm(c(V(:, e(k,1)), :) - c(V(:, e(k,2)), :), 2, 2);
end
q = 6*sqrt(2)*vol ./ max(h, [], 2).^3;
end
