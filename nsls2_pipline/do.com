#rm -rf *.xml
./mat.com

rm -rf 0*

mkdir 0.25
mkdir 0.30
mkdir 0.35
mkdir 0.40
mkdir 0.45
mkdir 0.50
mkdir 0.55
mkdir 0.60
mkdir 0.65
mkdir 0.70
mkdir 0.75

echo 'The "do.com" script simmulates what the database will do, and eventually should be deleted'

python do_sad.py << EOF
/GPFS/CENTRAL/XF17ID2/jjakoncic/Aug_15_2016/projID/XtalSamp_6_3/10/XtalSamp_6_3_00004.cbf
EOF
