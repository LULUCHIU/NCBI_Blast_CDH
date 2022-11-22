water -asequence query_tmp.fasta -bsequence target.fasta -gapopen 10 -gapextend 0.5 -aglobal3 Y -aglobal_outfile result
perl trim_m.pl result 
