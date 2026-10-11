% TEST_MG3D_HELMHOLTZ
% Multigrid for the CQ frequency-domain problem  -Lap u + (s/c)^2 u = f,
% s complex with Re s > 0, i.e. Helmholtz with complex wavenumber k = i s/c,
% homogeneous Neumann conditions, on the Gmsh box hierarchy.
%
% Sweep over |k| h_fine (fine-level resolution; k h ~ 0.6-1 is the 1-5 %
% accuracy range for P2) and arg(s) (tan(arg s) = Re k / Im k: small arg
% is a damped, nearly elliptic problem; arg -> 90 deg is the undamped
% Helmholtz limit). For each (|k|h, arg s) and each smoother ('sgs', 'ilu')
% the V-cycle contraction factor and the GMRES iteration count with the
% V-cycle as left preconditioner are recorded. Stiffness and mass matrices
% are assembled once per level; A(s) = S + (s/c)^2 M.
%
% Run from the repository root with addpath(genpath(pwd)).
%
% Last modified: 2026-10-11
% Changelog:
%   2026-10-11  First version.

clear; clc;

meshFile = fullfile('meshes', 'box_lc3_p2.msh');
nLevels  = 3;                  % 2 for a quick pass (~50k dofs), 3 for ~330k
c        = 1;
argDeg   = [10 45 70 85];      % arg(s) in degrees
khList   = [0.1 0.2 0.4 0.6 0.8 1.0 1.5 2.0];   % |k| h on the finest level
smoothers = {'sgs', 'ilu'};
tolGMRES = 1e-6;  maxitGMRES = 100;

obj0 = FEM3Dclass(meshFile);
levels = MG3D.buildHierarchy(obj0, nLevels);
L = nLevels;
S = cell(1, L); M = cell(1, L);
for k = 1:L
    S{k} = levels{k}.femStressMatrix();
    M{k} = levels{k}.femMassMatrix(@(x,y,z) 1 + 0*x);
end
mg = MG3D(levels, @(obj) obj.femStressMatrix(), 'Smoother', 'sgs', 'nu1', 2, 'nu2', 2);

% fine mesh size: median longest edge
Tf = levels{L}.mesh; cF = Tf.coord; V = Tf.ttrh(:, 1:4);
e = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; h = zeros(size(V,1), 6);
for j = 1:6, h(:, j) = vecnorm(cF(V(:,e(j,1)),:) - cF(V(:,e(j,2)),:), 2, 2); end
hFine = median(max(h, [], 2));
fprintf('%d levels, %d dofs on the finest level, h_fine = %.3f\n', L, size(S{L},1), hFine);

% right-hand side: smooth source
f  = @(x,y,z) exp(-((x-20).^2 + (y-3).^2 + (z-75).^2)/50);
b  = levels{L}.femSourceTerm(f);

rho  = nan(numel(khList), numel(argDeg), numel(smoothers));
its  = nan(size(rho));
flag = nan(size(rho));
tsol = nan(size(rho));
mark = {' ', '*'};
for ia = 1:numel(argDeg)
    for ik = 1:numel(khList)
        kabs = khList(ik)/hFine;                       % |k| = |s|/c
        s    = c*kabs*exp(1i*argDeg(ia)*pi/180);       % Re s > 0
        for k = 1:L
            mg.A{k} = S{k} + (s/c)^2 * M{k};
        end
        for is = 1:numel(smoothers)
            mg.Smoother = smoothers{is};
            mg.setup();
            r = mg.contractionFactor(12);
            [~, info] = mg.solve(b, 'Solver', 'gmres', 'Tol', tolGMRES, 'MaxIt', maxitGMRES);
            rho(ik, ia, is) = r;  its(ik, ia, is) = info.iter;
            flag(ik, ia, is) = info.flag;  tsol(ik, ia, is) = info.time;
            fprintf('arg %2d deg  |k|h %.2f  Re(k)/Im(k) %6.2f  %-4s  rho %7.3f  gmres %3d it%s  %.1f s\n', ...
                argDeg(ia), khList(ik), tan(argDeg(ia)*pi/180), smoothers{is}, r, info.iter, ...
                mark{1 + (info.flag ~= 0)}, info.time);
            
        end
    end
end

%% Summary table and plot
fprintf('\nGMRES iterations (tol %.0e); * = not converged in %d\n', tolGMRES, maxitGMRES);
for is = 1:numel(smoothers)
    fprintf('\n%s\n|k|h \\ arg(s): %s\n', smoothers{is}, sprintf('%8d', argDeg));
    for ik = 1:numel(khList)
        row = sprintf('%6.2f        ', khList(ik));
        for ia = 1:numel(argDeg)
            row = [row sprintf('%7d%s', its(ik, ia, is), ...
                  mark{1 + (flag(ik,ia,is) ~= 0)} )]; %#ok<AGROW>
        end
        fprintf('%s\n', row);
    end
end

figure;
for is = 1:numel(smoothers)
    subplot(1, numel(smoothers), is);
    semilogy(khList, squeeze(its(:, :, is)), 'o-'); grid on;
    xlabel('|k| h_{fine}'); ylabel('GMRES iterations');
    title(sprintf('%s, %d levels', smoothers{is}, L));
    legend(arrayfun(@(a) sprintf('arg s = %d^\\circ', a), argDeg, 'UniformOutput', false), 'Location', 'northwest');
end

if ~isfolder('results'), mkdir('results'); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmm'));
save(fullfile('results', sprintf('helmholtz_sweep_L%d_%s.mat', L, stamp)), ...
    'khList', 'argDeg', 'smoothers', 'rho', 'its', 'flag', 'tsol', 'hFine', 'nLevels', 'meshFile');
