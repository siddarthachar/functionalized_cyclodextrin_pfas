#!/bin/bash
shopt -s expand_aliases
source ~/.bashrc
prepgmx
export PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/bin:$PATH"
export LD_LIBRARY_PATH="/project2/andrewferguson/maxtopel/software/enhancements/plumed-2.8.1/installed-files/lib:$LD_LIBRARY_PATH"


# Input files
TPR="sys-prod.tpr"  # Replace with your TPR file
TRAJ="sys-prod.xtc"  # Replace with your trajectory file
OUT_DIR="processed_trajectory"

# Create output directory
mkdir -p $OUT_DIR

# Index group to isolate beta-CD
# Assuming beta-CD is named PRO in your topology
INDEX_FILE="bcd.ndx"
gmx make_ndx -f sys-neutral.gro -o temp.ndx << EOF
r PRO
q
EOF

grep -A 9999 "\[ PRO \]" temp.ndx > $INDEX_FILE # grep only under PRO

# Fit trajectory based on beta-CD atoms
# Select beta-CD group for fitting
# Select beta-CD group for output
gmx trjconv -s $TPR -f $TRAJ -o $OUT_DIR/fit_trajectory.xtc -fit rot+trans -n $INDEX_FILE << EOF
0  
0   
EOF

# Center trajectory to keep the system in the box
gmx trjconv -s $TPR -f $OUT_DIR/fit_trajectory.xtc -o $OUT_DIR/centered.xtc -pbc mol -center -n $INDEX_FILE << EOF
0   
0   
EOF

# Analyze RMSD relative to the starting structure
gmx rms -s $TPR -f $OUT_DIR/centered.xtc -o $OUT_DIR/rmsd.xvg -n $INDEX_FILE << EOF
0   
0   
EOF

# Analyze radius of gyration
# gmx gyrate -s $TPR -f $OUT_DIR/centered.xtc -o $OUT_DIR/gyrate.xvg -n $INDEX_FILE << EOF
# 0   
# EOF

# # Visualize hydrogen bonds
# gmx hbond -s $TPR -f $OUT_DIR/centered.xtc -num $OUT_DIR/hbonds.xvg -n $INDEX_FILE << EOF
# 0  
# 1   
# EOF

# # Analyze solvent accessible surface area (SASA)
# gmx sasa -s $TPR -f $OUT_DIR/centered.xtc -o $OUT_DIR/sasa.xvg -n $INDEX_FILE << EOF
# 0   # Select beta-CD group
# EOF

# Analyze distances between PFOS and beta-CD
# gmx distance -s $TPR -f $OUT_DIR/centered.xtc -oall $OUT_DIR/distances.xvg -n $INDEX_FILE << EOF
# 0   # Select beta-CD group
# 1   # Select PFOS group
# EOF

echo "Analysis complete. Results are in the $OUT_DIR directory."
