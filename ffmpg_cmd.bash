#!/usr/bin/env bash
# ======================================================================
# DESCRIPTION
#   Creates high quality movies that do not lose critical quality
#   after web-host re-encodes.
# 
#   Assumes 16:9 aspect ratio, but currently does not enforce --
#   implicit through bitrate
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
in_glob="fig_repo/prefix_%05d.png"
out_mov="new_movie.mp4"
height=1080

## Constants
IS_ODD=1

## Checks
(( fps & 1 == IS_ODD)) && {
  printf "fps must be even: $fps\n"  
  exit 150
} || :
gop=$(( fps / 2 ))
width=$(python -c "print( int(16/9.*$height) )")

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
  # --------------------------------------------------------------------
  # NOTE -- currently movie is not being forced to any size
  #      -- output will be the result of the size of the input frames.
  #      -- To change this, specify `size` or `pad` as a `-vf` flag.
  # --------------------------------------------------------------------
  # Padding to youtube 16:9 pixel size is recommended for frames that are
  # smaller--e.g. to 1920x1080p. See docs to change pad color.
# -vf pad=$width:$height:0:0
  -movflags +faststart
  "$out_mov"
)

ffmpeg "${opts[@]}"
