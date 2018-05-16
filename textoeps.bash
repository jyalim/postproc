#!/usr/bin/env bash
# ======================================================================
# textoeps -- script for compiling pgfplots/tikz into vector eps
# ----------------------------------------------------------------------
# Inputs:
#   in_tex  -- path to tex file (standalone)
#
# Returns:
#   None
#
# ======================================================================


tex_file="$1"
prefix=$(basename "${tex_file%.*}")

latex -halt-on-error -interaction=batchmode -jobname "${prefix}" "${tex_file}"
dvips -o "${prefix}.ps" "${prefix}.dvi"
ps2eps -f "${prefix}.ps"
rm "${prefix}.ps" "${prefix}.dvi" 
