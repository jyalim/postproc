#!/usr/bin/env bash
# ======================================================================
# DESCRIPTION
#   Creates high quality movies that do not lose critical quality
#   after web-host re-encodes.
# 
#   Assumes 16:9 aspect ratio, but currently does not enforce --
#   implicit through bitrate
# ----------------------------------------------------------------------
# USAGE
#     png_frames="path/to/frames/prefix_%04d.png"
#     mp4_movie="out.mp4"
#     framerate=60          
#     resize=600:600        
#     bash ffmpg_cmd.bash "$png_frames" "$mp4_movie" $framerate $resize
#  
#   * NOTE that the framerate is optional and defaults to 60 fps
#   * NOTE that resize is optional and defaults to the input frame size
#     and that resize must be in the format of 
#       width_pixels:height_pixels 
#     or to preserve aspect ratio:
#       width_pixels:-1 
#     which may not be as safe as -2 in place as -1 ( see the docs for 
#     more: https://trac.ffmpeg.org/wiki/Scaling )
#   * NOTE that resize is the 4th positional argument--framerate must be
#     specified if resize needs specification.
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
in_glob="${1:?Input files missing}"
out_mov="${2:?Output mp4 missing}"
fps="${3:-60}"
size="${4:-0}"
height=1080

## Constants
IS_ODD=1

## Checks
[[ $size != 0 ]] && {
  printf "Resize mode\nSetting size to: $size\n" 
  optional_vf_resize="-vf scale=$size"
  
}
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
  #      -- To change this, specify `scale` or `pad` as a `-vf` flag.
  # --------------------------------------------------------------------
  # Padding to youtube 16:9 pixel size is recommended for frames that are
  # smaller--e.g. to 1920x1080p. See docs to change pad color.
# -vf pad=$width:$height:0:0
  $optional_vf_resize
  -movflags +faststart
  "$out_mov"
)

ffmpeg "${opts[@]}"
