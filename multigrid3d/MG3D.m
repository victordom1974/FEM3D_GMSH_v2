classdef MG3D < handle
% MG3D  Geometric multigrid on a hierarchy of nested FEM3Dclass meshes.
%
%   mg = MG3D(levels, assembleFcn, Name, Value, ...)
%
%   levels      : cell array of FEM3Dclass objects, levels{1} the coarsest,
%                 levels{end} the finest, each obtained from the previous
%                 one with refineTtrh (see MG3D.buildHierarchy).
%   assembleFcn : function handle, A = assembleFcn(obj) returning the system
%                 matrix on the mesh of the FEM3Dclass object obj. Called on
%                 every level (rediscretization). With 'Galerkin', true the
%                 coarse matrices are P.'*A*P instead and assembleFcn is only
%                 called on the finest level.
%
%   Options (Name, Value):
%     'Smoother'  'sgs' (symmetric Gauss-Seidel, default), 'gs' (forward
%                 Gauss-Seidel), 'jacobi' (damped Jacobi), 'ilu' (ILU(0),
%                 robust to low-quality elements and local anisotropy)
%     'nu1','nu2' pre- and post-smoothing steps (default 2, 2)
%     'omega'     Jacobi damping (default 2/3)
%     'Galerkin'  false (default) | true
%     'Verbose'   false (default) | true
%
%   Methods:
%     x = mg.vcycle(b)            one V-cycle from x = 0 (as preconditioner)
%     x = mg.vcycle(b, x0)        one V-cycle from x0
%     [x, info] = mg.solve(b, ...) 'Solver': 'gmres' (default) | 'pcg' |
%                 'vcycle'; 'Tol' (1e-8), 'MaxIt' (200), 'Restart' ([] = no
%                 restart), 'x0'. info has flag, relres, iter, resvec, time.
%     rho = mg.contractionFactor(nIter)  asymptotic reduction factor of the
%                 V-cycle iteration on the homogeneous problem
%
%   Static:
%     [levels, P] = MG3D.buildHierarchy(obj0, nLevels)  refines obj0
%                 (FEM3Dclass) nLevels-1 times and builds the prolongations
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version: V-cycle, SGS/GS/Jacobi smoothers, direct
%               coarse solve, pcg/gmres interface, contraction factor.
%   2026-10-08  ILU(0) smoother.

    properties
        levels      % cell of FEM3Dclass, coarse -> fine
        A           % cell of system matrices per level
        P           % cell, P{k} prolongation from level k to level k+1
        nLevels
        Smoother = 'sgs'
        nu1 = 2
        nu2 = 2
        omega = 2/3
        Galerkin = false
        Verbose = false
    end
    properties (Access = private)
        Lw          % cell, tril(A{k}) for Gauss-Seidel
        Up          % cell, triu(A{k})
        Dinv        % cell, 1./diag(A{k}) for Jacobi
        Lilu        % cell, ILU(0) factors
        Uilu
        coarseFactor
    end

    methods
        function obj = MG3D(levels, assembleFcn, varargin)
            p = inputParser;
            addParameter(p, 'Smoother', 'sgs', @(s) any(strcmpi(s, {'sgs','gs','jacobi','ilu'})));
            addParameter(p, 'nu1', 2);
            addParameter(p, 'nu2', 2);
            addParameter(p, 'omega', 2/3);
            addParameter(p, 'Galerkin', false);
            addParameter(p, 'Verbose', false);
            parse(p, varargin{:});
            obj.Smoother = lower(p.Results.Smoother);
            obj.nu1 = p.Results.nu1;  obj.nu2 = p.Results.nu2;
            obj.omega = p.Results.omega;
            obj.Galerkin = p.Results.Galerkin;
            obj.Verbose = p.Results.Verbose;

            obj.levels  = levels(:).';
            obj.nLevels = numel(levels);
            L = obj.nLevels;

            % Prolongations
            obj.P = cell(1, L-1);
            for k = 1:L-1
                obj.P{k} = prolongation3D(levels{k}, levels{k+1}.mesh);
            end

            % Matrices
            obj.A = cell(1, L);
            if obj.Galerkin
                obj.A{L} = assembleFcn(levels{L});
                for k = L-1:-1:1
                    obj.A{k} = obj.P{k}.' * obj.A{k+1} * obj.P{k};
                end
            else
                for k = 1:L
                    obj.A{k} = assembleFcn(levels{k});
                end
            end
            obj.setup();
        end

        function setup(obj)
            % Smoother splittings and coarse factorization
            L = obj.nLevels;
            obj.Lw = cell(1, L); obj.Up = cell(1, L); obj.Dinv = cell(1, L);
            obj.Lilu = cell(1, L); obj.Uilu = cell(1, L);
            for k = 2:L
                switch obj.Smoother
                    case {'sgs', 'gs'}
                        obj.Lw{k} = tril(obj.A{k});
                        obj.Up{k} = triu(obj.A{k});
                    case 'jacobi'
                        obj.Dinv{k} = 1 ./ full(diag(obj.A{k}));
                    case 'ilu'
                        [obj.Lilu{k}, obj.Uilu{k}] = ilu(obj.A{k}, struct('type', 'nofill'));
                end
            end
            obj.coarseFactor = decomposition(obj.A{1});
            if obj.Verbose
                for k = 1:L
                    fprintf('MG3D level %d: %8d dofs, nnz %10d\n', k, size(obj.A{k},1), nnz(obj.A{k}));
                end
            end
        end

        function x = vcycle(obj, b, x0)
            if nargin < 3 || isempty(x0)
                x0 = zeros(size(b));
            end
            x = obj.cycle(obj.nLevels, b, x0);
        end

        function x = cycle(obj, k, b, x)
            if k == 1
                x = obj.coarseFactor \ b;
                return
            end
            x  = obj.smooth(k, x, b, obj.nu1, 'forward');
            r  = b - obj.A{k} * x;
            rc = obj.P{k-1}.' * r;
            ec = obj.cycle(k-1, rc, zeros(size(rc)));
            x  = x + obj.P{k-1} * ec;
            x  = obj.smooth(k, x, b, obj.nu2, 'backward');
        end

        function x = smooth(obj, k, x, b, nu, direction)
            Ak = obj.A{k};
            for it = 1:nu
                switch obj.Smoother
                    case 'jacobi'
                        x = x + obj.omega * (obj.Dinv{k} .* (b - Ak*x));
                    case 'gs'
                        x = x + obj.Lw{k} \ (b - Ak*x);
                    case 'ilu'
                        x = x + obj.Uilu{k} \ (obj.Lilu{k} \ (b - Ak*x));
                    case 'sgs'
                        if strcmp(direction, 'forward')
                            x = x + obj.Lw{k} \ (b - Ak*x);
                            x = x + obj.Up{k} \ (b - Ak*x);
                        else
                            x = x + obj.Up{k} \ (b - Ak*x);
                            x = x + obj.Lw{k} \ (b - Ak*x);
                        end
                end
            end
        end

        function [x, info] = solve(obj, b, varargin)
            p = inputParser;
            addParameter(p, 'Solver', 'gmres', @(s) any(strcmpi(s, {'gmres','pcg','vcycle'})));
            addParameter(p, 'Tol', 1e-8);
            addParameter(p, 'MaxIt', 200);
            addParameter(p, 'Restart', []);
            addParameter(p, 'x0', []);
            parse(p, varargin{:});
            tol = p.Results.Tol; maxit = p.Results.MaxIt; x0 = p.Results.x0;
            if isempty(x0), x0 = zeros(size(b)); end
            Af = obj.A{obj.nLevels};
            M  = @(r) obj.vcycle(r);
            t = tic;
            switch lower(p.Results.Solver)
                case 'gmres'
                    [x, flag, relres, iter, resvec] = gmres(Af, b, p.Results.Restart, tol, maxit, M, [], x0);
                    if isempty(p.Results.Restart), iter = iter(2); else, iter = (iter(1)-1)*p.Results.Restart + iter(2); end
                case 'pcg'
                    [x, flag, relres, iter, resvec] = pcg(Af, b, tol, maxit, M, [], x0);
                case 'vcycle'
                    x = x0; nb = norm(b);
                    resvec = zeros(maxit+1, 1); resvec(1) = norm(b - Af*x);
                    flag = 1;
                    for iter = 1:maxit
                        x = obj.vcycle(b, x);
                        resvec(iter+1) = norm(b - Af*x);
                        if resvec(iter+1) <= tol*nb, flag = 0; break; end
                    end
                    resvec = resvec(1:iter+1); relres = resvec(end)/nb;
            end
            info = struct('flag', flag, 'relres', relres, 'iter', iter, ...
                          'resvec', resvec, 'time', toc(t), 'solver', lower(p.Results.Solver));
            if obj.Verbose
                fprintf('MG3D %s: %d iterations, relres %.2e, %.2f s\n', info.solver, iter, relres, info.time);
            end
        end

        function [rho, hist] = contractionFactor(obj, nIter)
            % Asymptotic reduction factor of the V-cycle on A x = 0 from a
            % random start: ratio of consecutive error norms after transients.
            if nargin < 2, nIter = 15; end
            rng(0);
            n = size(obj.A{obj.nLevels}, 1);
            x = randn(n, 1);
            hist = zeros(nIter, 1);
            b = zeros(n, 1);
            for it = 1:nIter
                xn = obj.vcycle(b, x);
                hist(it) = norm(xn) / norm(x);
                x = xn;
            end
            rho = hist(end);
        end
    end

    methods (Static)
        function [levels, P] = buildHierarchy(obj0, nLevels)
            % Refine obj0 (FEM3Dclass) nLevels-1 times.
            levels = cell(1, nLevels);
            levels{1} = obj0;
            P = cell(1, nLevels-1);
            for k = 2:nLevels
                Tf = refineTtrh(levels{k-1}.mesh);
                levels{k} = FEM3Dclass(Tf);
                P{k-1} = prolongation3D(levels{k-1}, Tf);
            end
        end
    end
end
