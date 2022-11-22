mmseqs createdb query.fasta queryDB
mmseqs createdb target.fasta targetDB
mmseqs search queryDB targetDB resultDB tmp -s 11 --alignment-mode 0 --gap-open 15 --gap-extend 2 --min-aln-len 20 -a
mmseqs convertalis queryDB targetDB resultDB resultDB.m8
#mmseqs prefilter queryDB targetDB resultDB_pref -s 11
#mmseqs align queryDB targetDB resultDB_pref resultDB_aln
#mmseqs createtsv queryDB targetDB resultDB_aln resultDB_aln.txt
