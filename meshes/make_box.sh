#!/usr/bin/env bash
# Generate the test box mesh with Gmsh (needs the official binary with Netgen).
# Usage: meshes/make_box.sh [lc] [order]     (defaults: lc=3, order=2)
# Last modified: 2026-10-11
lc=${1:-3}; order=${2:-2}
gmsh -3 -order "$order" -format msh41 \
  -setnumber x0 10 -setnumber y0 1.5 -setnumber z0 0 \
  -setnumber Lx 20 -setnumber Ly 3 -setnumber Lz 150 \
  -setnumber lc "$lc" "$(dirname "$0")/box.geo" -o "$(dirname "$0")/box_lc${lc}_p${order}.msh"
