#This script should take a reference genome in fasta format and produce de novo repeat annotations. Outputs will include be identified repeat families, a full annotation of repeats and summary tables and graphs of repeat annotations  
##First input should the reference; second input should be the desired working directory
REFERENCE=$1
wDIR=PATH_TO_WORKING_DIRECTORY_FOR_ANALYSIS
REPEATMODELER2=PATH_TO_REPEATMODELER2_DIR
REPEATMASKER=PATH_TO_REPEATMASKER_DIR

rm -r $wDIR/RepeatModeler2
mkdir -p $wDIR/RepeatModeler2

#Make Database
cd $wDIR/RepeatModeler2
rnREF=`echo $REFERENCE | perl -pe 's#.*/##' | perl -pe 's#.gz$##'`
dbNAME=`echo $REFERENCE | perl -pe 's#.*/##' | perl -pe 's#\.f.*.gz$##'`
zcat $REFERENCE > $wDIR/RepeatModeler2/$rnREF
#Build Database
$REPEATMODELER2/BuildDatabase -engine ncbi -name $dbNAME $wDIR/RepeatModeler2/$rnREF

##RepeatModeler
$REPEATMODELER2/RepeatModeler -threads 24 -database $dbNAME -engine ncbi -genomeSampleSizeMax 243000000 -LTRStruct

rm -r $wDIR/RepeatModeler2/RM_*
echo 'RepeatModeler2 Complete'

rm -r $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7
mkdir -p $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7

cd $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7
cp $wDIR/RepeatModeler2/$dbNAME-families.fa $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$dbNAME-families.fa
cp $wDIR/RepeatModeler2/$rnREF $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF
$REPEATMASKER/RepeatMasker -e rmblast -gccalc -s -a -pa 24 -lib  $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$dbNAME-families.fa $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF

GENOMESIZE=`grep -v ">" $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF | perl -pe -chomp | wc -m`

perl $REPEATMASKER/util/calcDivergenceFromAlign.pl -s $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.divsum -a  $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.GC-Adjusted.align $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.align
perl $REPEATMASKER/util/createRepeatLandscape.pl -div $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$rnREF.divsum -t "$dbNAME Repeat Landscape" -g $GENOMESIZE > $wDIR/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/$dbNAME.repeat_landscape.html
echo 'RepeatMasker and Graphs Complete'
