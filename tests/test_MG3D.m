% TEST_MG3D
% Geometric multigrid for  -Lap u + u = f  on a box with homogeneous Neumann
% conditions. Exact solution: product of cosines with zero normal derivative
% on the six faces.
%
% Mesh: a box meshed by Gmsh through gmshBox (optimized, no slivers). If
% Gmsh is not available the test falls back to Examples/Elasticity/beam.msh,
% whose four slivers (quality 0.012) are inherited by every level and
% degrade the point smoothers (rho ~ 0.8 instead of ~0.1); the ILU(0)
% smoother is much less sensitive to this.
%
% Reports: mesh quality, discretization error, Galerkin vs rediscretization
% gap, V-cycle contraction factor and iteration counts for the smoothers
% 'sgs' and 'ilu', and the time of a direct solve for reference.
% Run from the repository root with addpath(genpath(pwd)).
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version.
%   2026-10-08  Gmsh box via gmshBox, sgs vs ilu comparison, quality report.

clear; clc;

deg     = 2;
nLevels = 3;
a = [10 1.5 0]; Lbox = [20 3 150];           % same box as beam.msh

try
    [obj0, ginfo] = gmshBox(Lbox, 3, deg, 'Origin', a, 'Verbose', true);
    meshName = 'gmshBox';
catch err
    fprintf('gmshBox not available (%s); using beam.msh\n', err.message);
    obj0 = FEM3Dclass(fullfile('Examples', 'Elasticity', 'beam.msh'));
    meshName = 'beam.msh';
end
q0 = meshQuality3D(obj0.mesh);
fprintf('coarse mesh %s: %d ttrh, %d nodes, quality min %.3f, median %.3f, #(q<0.1) = %d\n', ...
    meshName, size(obj0.mesh.ttrh,1), obj0.mesh.nNodes, min(q0), median(q0), nnz(q0 < 0.1));

uex = @(x,y,z) cos(pi*(x-a(1))/Lbox(1)) .* cos(pi*(y-a(2))/Lbox(2)) .* cos(pi*(z-a(3))/Lbox(3));
cf  = 1 + pi^2*sum(1./Lbox.^2);
f   = @(x,y,z) cf*uex(x,y,z);
assembleFcn = @(obj) obj.femStressMatrix() + obj.femMassMatrix(@(x,y,z) 1 + 0*x);

tic; levels = MG3D.buildHierarchy(obj0, nLevels); fprintf('hierarchy built in %.2f s\n', toc);

objF = levels{end};
b  = objF.femSourceTerm(f);
c  = objF.mesh.coord;
uN = uex(c(:,1), c(:,2), c(:,3));

results = struct('smoother', {}, 'rho', {}, 'pcg', {}, 'gmres', {}, 'tSetup', {}, 'tGmres', {});
for sm = {'sgs', 'ilu'}
    tic;
    mg = MG3D(levels, assembleFcn, 'Smoother', sm{1}, 'nu1', 2, 'nu2', 2);
    tSetup = toc;
    if strcmp(sm{1}, 'sgs')
        Af = mg.A{end};
        G  = mg.P{end}.' * Af * mg.P{end};
        fprintf('Galerkin vs rediscretization (level %d): |P''AP - A_c|/|A_c| = %.2e\n', ...
            nLevels-1, norm(G - mg.A{end-1}, 1)/norm(mg.A{end-1}, 1));
        tic; ud = Af \ b; tDirect = toc;
        fprintf('direct solve: %d dofs, %.2f s, nodal error vs exact %.2e\n', ...
            numel(b), tDirect, norm(ud - uN, inf)/norm(uN, inf));
    end
    rho = mg.contractionFactor(15);
    [~, ip] = mg.solve(b, 'Solver', 'pcg',   'Tol', 1e-10, 'MaxIt', 100);
    [u, ig] = mg.solve(b, 'Solver', 'gmres', 'Tol', 1e-10, 'MaxIt', 100);
    fprintf('%-4s: setup %.2f s, rho %.3f, pcg %3d it (%.2f s), gmres %3d it (%.2f s), |u - u_direct| %.1e\n', ...
        sm{1}, tSetup, rho, ip.iter, ip.time, ig.iter, ig.time, norm(u - ud, inf)/norm(ud, inf));
    results(end+1) = struct('smoother', sm{1}, 'rho', rho, 'pcg', ip.iter, ...
        'gmres', ig.iter, 'tSetup', tSetup, 'tGmres', ig.time); %#ok<SAGROW>
end
