mmseqs createdb query.fasta queryDB
mmseqs createdb target.fasta targetDB
mmseqs prefilter queryDB targetDB resultDB_pref -s 11
mmseqs align queryDB targetDB resultDB_pref resultDB_aln
mmseqs createtsv queryDB targetDB resultDB_aln resultDB_aln.txt
