makeblastdb -in query.fasta -out queryDB -dbtype prot
blastp -query target.fasta -db queryDB -outfmt 7 -out result.txt
grep -v "#" result.txt > result_final.txt
