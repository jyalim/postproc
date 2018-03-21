#!/usr/bin/env bash
# ======================================================================
# gsgs -- ghostscript grayscale
# ----------------------------------------------------------------------
# Converts color pdf to grayscale
# ----------------------------------------------------------------------
# Inputs:
#   in_pdf  -- path to color pdf
#   out_pdf -- OPTIONAL (default appends _bw to pdf prefix)
#
# Returns:
#   None
#
# From:
#     https://superuser.com/questions/104656/convert-a-pdf-to-greyscale-on-the-command-line-in-floss
# ======================================================================

in_pdf="${1}"
out_pdf="${2:-${in_pdf/.pdf/_bw.pdf}}"

opts=(
  -sOutputFile="$out_pdf"
  -sDEVICE=pdfwrite  
  -sColorConversionStrategy=Gray  
  -dProcessColorModel=/DeviceGray  
  -dCompatibilityLevel=1.4  
  -dNOPAUSE  
  -dAutoRotatePages=/None
  -dBATCH  
)

gs ${opts[@]} "$in_pdf"

