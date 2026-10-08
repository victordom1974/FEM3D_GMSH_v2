function [obj, info] = gmshBox(dims, lc, deg, varargin)
% GMSHBOX  Mesh a box with Gmsh (meshes/box.geo) and load it as FEM3Dclass.
%
%   [obj, info] = gmshBox(dims, lc, deg, Name, Value, ...)
%
%   dims : [Lx Ly Lz] box dimensions
%   lc   : target element size
%   deg  : polynomial degree of the Lagrange elements (1..4)
%
%   Options:
%     'Origin'  [x0 y0 z0], default [0 0 0]
%     'Output'  .msh file name, default fullfile(tempdir, 'box.msh')
%     'Gmsh'    gmsh executable, default 'gmsh'
%     'Verbose' false (default) | true (print the Gmsh log)
%
%   info.cmd      command line used
%   info.quality  per-tetrahedron quality (see meshQuality3D)
%   info.log      Gmsh output
%
%   The .geo switches on mesh optimization (Netgen), so the minimum
%   quality should be well above 0.2. Check info.quality: the quality of
%   this mesh bounds the quality of every refined level.
%
% Last modified: 2026-10-08
% Changelog:
%   2026-10-08  First version.

p = inputParser;
addParameter(p, 'Origin', [0 0 0]);
addParameter(p, 'Output', fullfile(tempdir, 'box.msh'));
addParameter(p, 'Gmsh', 'gmsh');
addParameter(p, 'Verbose', false);
parse(p, varargin{:});
o = p.Results;

geo = fullfile(fileparts(mfilename('fullpath')), '..', 'meshes', 'box.geo');
if ~isfile(geo)
    error('gmshBox:geo', 'box.geo not found at %s', geo);
end

[st, ver] = system(sprintf('"%s" --version', o.Gmsh));
if st ~= 0
    error('gmshBox:nogmsh', 'Gmsh executable "%s" not found (%s).', o.Gmsh, strtrim(ver));
end

cmd = sprintf(['"%s" -3 -order %d -format msh41 ' ...
               '-setnumber x0 %.15g -setnumber y0 %.15g -setnumber z0 %.15g ' ...
               '-setnumber Lx %.15g -setnumber Ly %.15g -setnumber Lz %.15g ' ...
               '-setnumber lc %.15g "%s" -o "%s" -v 2'], ...
               o.Gmsh, deg, o.Origin(1), o.Origin(2), o.Origin(3), ...
               dims(1), dims(2), dims(3), lc, geo, o.Output);
[st, log] = system(cmd);
if o.Verbose, disp(log); end
if st ~= 0 || ~isfile(o.Output)
    error('gmshBox:failed', 'Gmsh failed:\n%s\n%s', cmd, log);
end

obj = FEM3Dclass(o.Output);
info.cmd     = cmd;
info.log     = log;
info.version = strtrim(ver);
info.quality = meshQuality3D(obj.mesh);
if o.Verbose
    fprintf('gmshBox: Gmsh %s, %d tetrahedra, %d nodes, quality min %.3f median %.3f\n', ...
        info.version, size(obj.mesh.ttrh,1), obj.mesh.nNodes, min(info.quality), median(info.quality));
end
end
