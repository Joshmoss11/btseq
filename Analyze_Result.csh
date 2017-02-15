#!/bin/tcsh -f

if ($#argv != 2) then
	echo "Usage: $0 <fastq file name> <Fasta dir>"
	exit 0
endif

set FastqName = $1
set FastaDir = $2

echo "Running bismark_methylation_extractor"
set sam_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print $NF}'`_bismark_bt2.sam
/cs/icore/joshua.moss/scripts/bismark/bismark_v0.13.1/bismark_methylation_extractor -s $sam_name --bedGraph --counts --cytosine_report --genome_folder ${FastaDir} --CX --comprehensive -o ${FastaDir}/Output

# Format the summary report nicely for all C's
set CX_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print $NF}'`_bismark_bt2.CX_report.txt

echo "Generating Summary and Hist files"

set report_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print $NF}'`_bismark_bt2_SE_report.txt
cat $report_name | grep "Sequences analysed in total" | awk '{print "Total # of reads:\t" $NF}' > $CX_name.summary
echo "Total reads aligned:\t`cat $sam_name | grep -v -P '^@' | wc -l`" >> $CX_name.summary

foreach tname (`cat ${FastaDir}/*.fa | grep ">" | sed s/">"// | awk '{print $1}'`)
echo "Read aligned to ${tname}:\t`cat $sam_name | grep -v -P '^@' | grep ${tname} | wc -l`" >> $CX_name.summary
end

foreach tname (`cat ${FastaDir}/*.fa | grep ">" | sed s/">"// | awk '{print $1}'`)
echo >> $CX_name.summary
echo $tname >> $CX_name.summary
echo "CpG sites" >> $CX_name.summary
echo "Position\t#C\t#T\t%Meth" >> $CX_name.summary
cat $CX_name | grep ${tname} | grep + | awk 'BEGIN {OFS="\t"; print "Position","Context","C_count","T_count","MethLevel"} {print $2,substr($7,0,2),$4,$5,($5>0)?100*$4/($4+$5):0}' > $CX_name.summary_all
cat $CX_name.summary_all | grep -w CG | awk '{OFS="\t"; print $1,$3,$4,$5}' >> $CX_name.summary
echo "\nCpA sites" >> $CX_name.summary
echo "Position\t#C\t#T\t%Meth" >> $CX_name.summary
cat $CX_name.summary_all | grep -w CA | awk '{OFS="\t"; print $1,$3,$4,$5}' >> $CX_name.summary
echo "\nAll C sites" >> $CX_name.summary
cat $CX_name.summary_all >> $CX_name.summary
rm $CX_name.summary_all

# Output the pattern distribution CpG
set CpGcontext_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print "CpG_context_" $NF}'`_bismark_bt2.txt
set CpGList = `cat $CX_name | grep $tname | grep + | grep -w CG | awk '{print $2}' | tr '\n' '\t' | sed 's/\t$//'`
cat $CpGcontext_name | grep ${tname} | awk '{FSO=" ";print $1, $4 ":" $2}' | awk '{ stuff[$1] = stuff[$1] $2 " " } END { for( s in stuff ) print s, stuff[s];}' | cut -d " " -f 2- | sed 's/^ *//' | tr ' ' ':' | awk -F':' -v list="$CpGList" 'BEGIN {split(list,arr," ");} {for (i=1;i<=length(arr);i++) b[arr[i]]="-"; for (i=2; i<=NF; i+=2) $(i)=="+"?b[$(i-1)]="C":b[$(i-1)]="T"; for (i=1;i<=length(arr);i++) printf b[arr[i]] "\t"; printf "\n"}' | sort | uniq -c | sort -nr | awk -v list="$CpGList" 'BEGIN {printf "Count\t" list "\n"} {print $0}' | sed 's/^ *//' | tr ' ' '\t' > $CpGcontext_name.hist.tmp
touch $CpGcontext_name.hist
echo $tname >> $CpGcontext_name.hist
cat $CpGcontext_name.hist.tmp | awk 'BEGIN {A=0;C=0;G=0;T=0;OFS="\t"} {if (NR==1) {printf $1 "\t" "# of A\t# of C\t# of G\t# of T\t# of sites\t"; for (i=2; i<=NF; i++) {printf $i "\t"} printf "\n"; next} for (i=2; i<=NF; i++) {if ($i=="A") {A=A+1} if ($i=="C") {C=C+1} if ($i=="G") {G=G+1} if ($i=="T") {T=T+1}} printf $1 "\t" A "\t" C "\t" G "\t" T "\t" (A+C+G+T) "\t"; for (i=2; i<=NF; i++) {printf $i "\t"} printf "\n";A=0;C=0;G=0;T=0}' >> $CpGcontext_name.hist
echo >> $CpGcontext_name.hist
rm $CpGcontext_name.hist.tmp

# Output the pattern distribution CpA
set CpAcontext_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print "CpA_context_" $NF}'`_bismark_bt2.txt
set CHHcontext_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print "CHH_context_" $NF}'`_bismark_bt2.txt
set CHGcontext_name = ${FastaDir}/Output/`echo ${FastqName} | awk -F"/" '{print "CHG_context_" $NF}'`_bismark_bt2.txt
set CpAList = `cat $CX_name | grep $tname | grep + |  grep -P 'CA.$' | awk '{print $2}' | tr '\n' '\t' | sed 's/\t$//'`
cat ${CHHcontext_name} ${CHGcontext_name} | grep ${tname} | awk '{FSO=" ";print $1, $4 ":" $2}' | awk '{ stuff[$1] = stuff[$1] $2 " " } END { for( s in stuff ) print s, stuff[s];}' | cut -d " " -f 2- | sed 's/^ *//' | tr ' ' ':' | awk -F':' -v list="$CpAList" 'BEGIN {split(list,arr," ");} {for (i=1;i<=length(arr);i++) b[arr[i]]="-"; for (i=2; i<=NF; i+=2) $(i)=="+"?b[$(i-1)]="C":b[$(i-1)]="T"; for (i=1;i<=length(arr);i++) printf b[arr[i]] "\t"; printf "\n"}' | sort | uniq -c | sort -nr | awk -v list="$CpAList" 'BEGIN {printf "Count\t" list "\n"} {print $0}' | sed 's/^ *//' | tr ' ' '\t' > $CpAcontext_name.hist.tmp
touch $CpAcontext_name.hist
echo $tname >> $CpAcontext_name.hist
cat $CpAcontext_name.hist.tmp | awk 'BEGIN {A=0;C=0;G=0;T=0;OFS="\t"} {if (NR==1) {printf $1 "\t" "# of A\t# of C\t# of G\t# of T\t# of sites\t"; for (i=2; i<=NF; i++) {printf $i "\t"} printf "\n"; next} for (i=2; i<=NF; i++) {if ($i=="A") {A=A+1} if ($i=="C") {C=C+1} if ($i=="G") {G=G+1} if ($i=="T") {T=T+1}} printf $1 "\t" A "\t" C "\t" G "\t" T "\t" (A+C+G+T) "\t"; for (i=2; i<=NF; i++) {printf $i "\t"} printf "\n";A=0;C=0;G=0;T=0}' >> $CpAcontext_name.hist
echo >> $CpAcontext_name.hist
rm $CpAcontext_name.hist.tmp



end




