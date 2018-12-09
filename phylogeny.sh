
cat phage.fasta <(awk -F\| '/>/{$0=">"$2}1' uniprot.fasta) orfs/orf-nbla.fa > nbla.fasta

samtools faidx nbla.fasta
cut -f1 nbla.fasta.fai | grep -vf outliers | xargs samtools faidx nbla.fasta | MAFFT_BINARIES= mafft --localpair --maxiterate 1000 - > nbla.mafft.fasta
# MAFFT_BINARIES= mafft --localpair --maxiterate 1000 nbla.fasta > nbla.mafft.fasta

# /opt/EvalMSA/bin/EvalMSA nbla.woutliers.mafft.fasta /opt/EvalMSA/res/Matrix/blosum62
# gawk '/Median/{med=$2}!match($0, / (.*)\(([0-9]+)\)/, a){next} $0~/gene/ && a[2] > thr { thr = a[2] } $0!~/gene/ && (a[2] > thr || a[2] < med) { print a[1],a[2] }' nbla_results.txt > uniprot.outliers
# gawk 'NR==FNR{o[">"$1]=1;next} />/{p=1} o[$1] {p=0} p' uniprot.outliers nbla.woutliers.fasta > nbla.fasta
# mafft nbla.fasta > nbla.mafft.fasta

mkdir -p phy
trimal -in nbla.mafft.fasta -phylip -automated1 > phy/nbla.phy

cat << PF > phy/partition_finder.cfg
alignment = nbla.phy;
branchlengths = linked;
models = all;
model_selection = aicc;
[data_blocks]
nblA = 1-$(awk 'NR>1{exit}{print$2}' phy/nbla.phy);
[schemes]
search = greedy;
PF

/opt/partitionfinder/PartitionFinderProtein.py phy -p 10 -q --raxml --rcluster-max=100
model=$(awk '$1==1{print tolower($3)}' phy/analysis/best_scheme.txt | cut -d+ -f1)
OMP_NUM_THREADS=10 FastTreeMP "-$model" -gamma phy/nbla.phy > phy/nbla.nwk
