#!/usr/bin/env bash
# ======================================================================
# epstocpdf -- eps to cropped pdf
# ----------------------------------------------------------------------
# Converts eps to cropped pdf
# ----------------------------------------------------------------------
# Inputs:
#   in_eps  -- path to input eps, MUST HAVE .eps extension
#
# Returns:
#   None
#
# ======================================================================

in_eps="${1}"
bn="${1#*.eps}"
tmp_pdf="${bn}-tmp.pdf"
out_pdf="${bn}.pdf"

epstopdf "$in_eps" "$tmp_pdf"   && \
  pdfcrop "$tmp_pdf" "$out_pdf" && \
  rm "$tmp_pdf"                 || {
  echo 'FAILURE' 
}


