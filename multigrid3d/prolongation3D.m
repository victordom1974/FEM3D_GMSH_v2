function P = prolongation3D(coarse, Tf)
% PROLONGATION3D  Nested prolongation between consecutive refinement levels.
%
%   P = prolongation3D(coarse, Tf)
%
%   coarse : FEM3Dclass object (or mesh struct) of the coarse level
%   Tf     : mesh struct of the fine level, produced by refineTtrh(coarse),
%            so that Tf.refinement holds nodeParent and nodeBaryc
%
%   P      : sparse (Tf.nNodes x coarse.nNodes). P*uc is the fine nodal
%            vector of the coarse FE function with nodal values uc; since the
%            meshes are nested, P*uc represents exactly the same function.
%            Restriction: P.' (Galerkin). Coarse operator: P.'*A*P, or a
%            rediscretization on the coarse mesh (identical for nested
%            meshes and exact integration).
%
%   The coarse local basis Nj3D is evaluated at the barycentric coordinates
%   of each fine node inside its parent tetrahedron; no point location is
%   needed. Works for any degree 1-4 (both levels have the same degree).
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version (extracted from the tests).

if isa(coarse, 'FEM3Dclass')
    Tc  = coarse.mesh;
    Nj  = coarse.Nj3D;
else
    Tc  = coarse;
    [~, Nj] = FEM3Dclass.basisNj(Tc.deg);
end
if ~isfield(Tf, 'refinement') || isempty(Tf.refinement.nodeParent)
    error('prolongation3D:noHierarchy', 'Tf must come from refineTtrh (missing Tf.refinement).');
end
if size(Tf.refinement.nodeBaryc, 1) ~= Tf.nNodes || max(Tf.refinement.nodeParent) > size(Tc.ttrh, 1)
    error('prolongation3D:mismatch', 'Tf does not look like a refinement of the given coarse mesh.');
end

nL   = size(Tc.ttrh, 2);
nF   = Tf.nNodes;
b    = Tf.refinement.nodeBaryc;
vals = zeros(nF, nL);
for j = 1:nL
    vals(:, j) = Nj{j}(b(:,1), b(:,2), b(:,3), b(:,4));
end
cols = Tc.ttrh(Tf.refinement.nodeParent, :);
rows = repmat((1:nF).', 1, nL);
keep = abs(vals) > 1e-14;
P = sparse(rows(keep), cols(keep), vals(keep), nF, Tc.nNodes);
end
