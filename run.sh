
cat genes.fasta uniprot.fasta > nbla.woutliers.fasta
mafft nbla.woutliers.fasta > nbla.woutliers.mafft.fasta
/opt/EvalMSA/bin/EvalMSA nbla.woutliers.mafft.fasta /opt/EvalMSA/res/Matrix/blosum62
awk '/Median/{med=$2}!match($0, / (.*)\(([0-9]+)\)/, a){next} $0~/gene/ && a[2] > thr { thr = a[2] } $0!~/gene/ && (a[2] > thr || a[2] < med) { print a[1],a[2] }' nbla_results.txt
awk 'NR==FNR{o[">"$1]=1;next} />/{p=1} o[$1] {p=0} p' uniprot.outliers nbla.woutliers.fasta > nbla.fasta
mafft nbla.fasta > nbla.mafft.fasta

mkdir -p phy
trimal -in nbla.mafft.fasta -out phy/nbla.phy -phylip -automated1

cat <<- PF > phy/partition_finder.cfg
	# ALIGNMENT FILE #
	alignment = nbla.phy;

	# BRANCHLENGTHS #
	branchlengths = linked;

	# MODELS OF EVOLUTION #
	models = all;
	model_selection = aicc;

	# DATA BLOCKS #
	[data_blocks]
	nblA = 1-$(awk 'NR>1{exit}{print$2}' phy/nbla.phy);

	# SCHEMES #
	[schemes]
	search = greedy;
PF
