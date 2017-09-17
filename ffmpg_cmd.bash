#!/usr/bin/env bash
# ======================================================================
# ASSUMES 
#   16:9 aspect ratio
# ----------------------------------------------------------------------
# DESCRIPTION
#   Creates high quality movies that will not lose quality when
#   web-host re-encodes.
# ----------------------------------------------------------------------
# NOTE
#   After publishing to web-host (e.g. youtube), it may take some time
#   before high-quality stream is available. 
# ----------------------------------------------------------------------
# RESOURCES
# [1] - youtube's recommended settings
#   https://support.google.com/youtube/answer/1722171?hl=en-GB
# [2] - user `c0nw0nk` trying to determine ffmpeg flags for youtube
#   https://ffmpeg.zeranoe.com/forum/viewtopic.php?t=2318
# ======================================================================

## INPUTS
fps=100
in_glob="fig_aps_mov/B13000e-4_N2e4_Pr1e0_F30e-2_m256_tr5e3_full_ps_%05d_T.png"
out_mov="B13t.mp4"
height=1080
width=$(python -c "print( int(16/9.*$height) )")

## Constants
IS_ODD=1

(( fps & 1 == IS_ODD)) && {
  printf "fps must be even: ${fps}\n"  
  exit 150
} || :
gop=$((100/2))

## Subroutines
get_bit_opts() {
  # Rates are based on [1], youtube's recommended bitrates 
  python << __EOF 
height = '$height'
lookup = {
  '360' : {'min': '1000k','max': '1500k','buf':  '3000k'},
  '480' : {'min': '2500k','max': '4000k','buf':  '8000k'},
  '720' : {'min': '5000k','max': '7500k','buf': '15000k'},
  '1080': {'min': '8000k','max':'12000k','buf': '24000k'},
  '1440': {'min':'16000k','max':'24000k','buf': '48000k'},
  '2160': {'min':'35000k','max':'68000k','buf':'136000k'},
}
out_str = '-b:v {max:s} -minrate {min:s} -maxrate {max:s} -bufsize {buf:s}'
print( out_str.format(**lookup[height]) )
__EOF
}

opts=(
  -framerate $fps
  -i "$in_glob"
	-c:v libx264
	-flags:v "+cgop" 
  -g $gop
  -profile:v high -level 4.0
  $(get_bit_opts)
  -bf 2 
  -coder 1
  -crf 16
  -pix_fmt yuv420p
  -strict -2
  -r $fps
  # https://video.stackexchange.com/questions/9947/how-do-i-change-frame-size-preserving-width-using-ffmpeg
  -vf setsar=1:1 
  # padding to youtube 16:9 pixel size is recommended for frames that are
  # smaller--e.g. to 1920x1080p. See docs to change pad color.
# -vf pad=$width:$height:0:0
  -movflags +faststart
  "$out_mov"
)

ffmpeg "${opts[@]}"
