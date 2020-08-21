from Bio import SeqIO
from Bio.Seq import Seq
import os
import subprocess
from fuzzysearch import find_near_matches
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#Insert paths and file names for four lines below - lines 11,12,13,15

output_file_name = "Insert_path_with_file"
subprocess.call("makeblastdb -in " + "Insert_blast_db_fasta_promoters" + " -out ***Insert_output_file***" + " -dbtype nucl", shell=True)
subprocess.call("makeblastdb -in " + "Insert_blast_db_fasta_Transposases" + " -out ***Insert_output_file***" + " -dbtype nucl", shell=True)
#The block below looks for the part6 sequence before the barcode in the query fasta sequences and gets the barcode sequence. Also annotates if it finds the part6 sequence or not.
with open("Insert_query_fasta_file", "r") as input_handle:
	for record in SeqIO.parse(input_handle, "fasta"):
		my_s = record.seq
		bar_match=list(find_near_matches('TAATACGACTCACTATAGGGGATAGATGTCCACGAGCTCTCT', my_s, max_l_dist=4))
		if (len(bar_match)==0):
			barcode="no part6"
			part6="no part6"
		elif (len(bar_match)>1):
			barcode="more than 1 part6"
			part6="more than 1 part6"
		else:
			bar_start=bar_match[0].start
			bar_end=bar_match[0].end
			part6=my_s[bar_start:bar_end+40]
			barcode=my_s[bar_end:bar_end+20]
		print(part6)
		print(barcode)

#Insert paths and file names for the lines below - line35,37,38,41,42
#The block below looks at each fasta sequence opened and annotated for part6 above, and annotates which transposase and promoter sequence it contains. Finally outputs a table with all the annotations per read.
		with open("temp_file_name", 'w') as temp:
			SeqIO.write(record, temp, 'fasta')
		subprocess.call("blastn -num_threads 48 -query " + "temp_file_name" + " -db " + "tpase_blast_db_output_from_line13" + " -outfmt 5 -out " + "Tpase_results.xml", shell = True)
		Tpase_blast_records = NCBIXML.parse(open("Tpase_results.xml"))
		Tpase_blast_record = next(Tpase_blast_records)
		for Tpase_alignment in Tpase_blast_record.alignments:
			subprocess.call("blastn -num_threads 48 -query " + "temp_file_name" + " -db " + "promoters_blast_db_output_from_line12" + " -outfmt 5 -out " + "Promoter_results.xml", shell = True)
			blast_records = NCBIXML.parse(open("Promoter_results.xml"))
			blast_record = next(blast_records)
			for alignment in blast_record.alignments:
				with open(output_file_name, "a") as ofile:
					ofile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(record.description, alignment.hit_def, Tpase_alignment.hit_def, part6, barcode))
