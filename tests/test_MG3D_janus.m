% TEST_MG3D_JANUS_CQ
% CQ-type frequency problems on the Janus geometry (meshes/janus.geo):
%     -Lap u + (s^2 / c(x)^2) u = f / c(x)^2,   homogeneous Neumann,
% with c piecewise constant per physical volume (outer: 1; particle halves
% 'a': cA, 'b': cB) and s = sigma + i*omega, sigma > 0 fixed, omega swept.
% S (Laplacian) and Mc (mass with 1/c^2) are assembled once per level;
% A(s) = S + s^2 Mc. For each s: V-cycle contraction factor and GMRES
% iterations with the V-cycle (SGS 2+2) as preconditioner.
% Run from the repository root with addpath(genpath(pwd)). The hierarchy is
% reused between runs when meshFile and nLevels are unchanged.
%
% Last modified: 2026-10-12
% Changelog:
%   2026-10-12  First version.

clc;

meshFile = fullfile('meshes', 'janus_p1.msh');
nLevels  = 3;
cA = 4;  cB = 0.25;                         % velocities in the particle halves
sigma    = 2.2;                             % Re s (CQ contour)
omegaOverSigma = [0 0.5 1 2 4 8 16 32];     % Im s / Re s; arg s = atan(.)
tolGMRES = 1e-6;  maxitGMRES = 100;

% --- hierarchy (cached)
rebuild = ~exist('levels', 'var') || ~exist('cacheMesh', 'var') || ...
    ~strcmp(cacheMesh, meshFile) || cacheLevels ~= nLevels;
if rebuild
    obj0 = FEM3Dclass(meshFile);
    tic; levels = MG3D.buildHierarchy(obj0, nLevels); fprintf('hierarchy built in %.2f s\n', toc);
    cacheMesh = meshFile; cacheLevels = nLevels;
end
L = nLevels; objF = levels{L};
for k = 1:L
    fprintf('  level %d: %8d ttrh, %8d dofs\n', k, size(levels{k}.mesh.ttrh,1), levels{k}.mesh.nNodes);
end

% --- coefficient by domain tag (outer = 1; j+1 -> cA; j+2 -> cB)
cOf  = @(d) 1 + (cA - 1)*(mod(d, 100) == 1 & d > 1) + (cB - 1)*(mod(d, 100) == 2);
invc2 = @(x,y,z,d) 1 ./ cOf(d).^2 + 0*x;

% --- matrices once per level
S = cell(1, L); Mc = cell(1, L);
tic;
for k = 1:L
    S{k}  = levels{k}.femStressMatrix();
    Mc{k} = levels{k}.femMassMatrix(invc2);
end
fprintf('S and Mc assembled on %d levels in %.1f s\n', L, toc);
mg = MG3D(levels, @(obj) obj.femStressMatrix(), 'Smoother', 'sgs', 'nu1', 2, 'nu2', 2);

% --- right-hand side f / c^2, smooth source centred in the box
f  = @(x,y,z,d) exp(-(x.^2 + y.^2 + z.^2)/2) ./ cOf(d).^2;
b  = objF.femSourceTerm(f);

% fine mesh size (median longest edge), for the |k| h column (c = 1 region)
Tf = objF.mesh; cF = Tf.coord; V = Tf.ttrh(:, 1:4);
e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; h = zeros(size(V,1), 6);
for j = 1:6, h(:, j) = vecnorm(cF(V(:,e(j,1)),:) - cF(V(:,e(j,2)),:), 2, 2); end
hFine = median(max(h, [], 2));

mark = {' ', '*'};
res = nan(numel(omegaOverSigma), 4);
fprintf('\nsigma = %.2f, h_fine = %.3f, c = 1 / %g / %g\n', sigma, hFine, cA, cB);
fprintf('%8s %8s %8s %8s %7s %6s %7s\n', 'omega', '|s|', 'arg(deg)', '|s|h', 'rho', 'gmres', 'time');
for io = 1:numel(omegaOverSigma)
    s = sigma*(1 + 1i*omegaOverSigma(io));
    for k = 1:L
        mg.A{k} = S{k} + s^2 * Mc{k};
    end
    mg.setup();
    rho = mg.contractionFactor(12);
    [~, info] = mg.solve(b, 'Solver', 'gmres', 'Tol', tolGMRES, 'MaxIt', maxitGMRES);
    res(io, :) = [imag(s), rho, info.iter, info.time];
    fprintf('%8.2f %8.2f %8.1f %8.2f %7.3f %5d%s %6.1fs\n', imag(s), abs(s), angle(s)*180/pi, ...
        abs(s)*hFine, rho, info.iter, mark{1 + (info.flag ~= 0)}, info.time);
end