shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH"
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"

rm -r segment_r
rm -r segment_z
rm -r segment_angle

mkdir segment_r
mkdir segment_z
mkdir segment_angle

cd segment_r
cp ../reweight/hills_r.out .
plumed sum_hills --hills hills_r.out --stride 100 --mintozero

cd ../segment_z
cp ../reweight/hills_z.out .
plumed sum_hills --hills hills_z.out --stride 100 --mintozero # change these later to hills_z.out

cd ../segment_angle
cp ../reweight/hills_angle.out .
plumed sum_hills --hills hills_angle.out --stride 100 --mintozero
