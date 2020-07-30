
source params.cfg

mkdir -p phy
awk -F\| '/>/{$0=">"$2}1' uniprot.fasta | cat phage.fasta - orfs/orf-nbla.fa | seqkit grep -vf outliers.txt | mafft --localpair --maxiterate 1000 - > phy/nbla.mafft.fasta

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

PartitionFinderProtein.py phy -p "$CPUs" -q --raxml --rcluster-max=100
model=$(awk '$1==1{print tolower($3)}' phy/analysis/best_scheme.txt | cut -d+ -f1)
OMP_NUM_THREADS="$CPUs" FastTreeMP "-$model" -gamma phy/nbla.phy > phy/nbla.nwk
