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


in_tex="$1"
prefix=$(basename "${in_tex%.*}")

latex -halt-on-error -interaction=batchmode -jobname "${prefix}" "${in_tex}"
dvips -o "${prefix}.ps" "${prefix}.dvi"
ps2eps -f "${prefix}.ps"
rm "${prefix}.ps" "${prefix}.dvi" 
