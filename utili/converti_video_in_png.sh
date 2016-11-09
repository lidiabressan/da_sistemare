_filee=$1
rad_file=${_filee//'.MPG'/''}
rad_file=${rad_file//'.mpg'/''}
rad_file=${rad_file//'MOV'/'mov'}
mkdir $rad_file
ffmpeg -i $_filee -f image2 "$rad_file/$rad_file"%05d.png
mv $_filee $rad_file/
