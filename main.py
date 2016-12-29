from Bio import SeqIO
import json
from itertools import count

def shred(seqrecord, k=33):
	kmer_set = set()
	for i in range(0,len(seqrecord.seq) - k):
		kmer_set.add(str(seqrecord.seq[i:i+k]))
	return seqrecord.name, kmer_set

def build_database(seqs, fp):
	db = dict(kmers=dict(), names=dict())
	for seq, i in zip(seqs, count()):
		name, kmer_set = shred(seq)
		#db[name] = kmer_set
		for kmer in kmer_set:
			if kmer not in db['kmers']:
				#sets can't be JSON serialized
				db['kmers'][kmer] = list()
			db['kmers'][kmer].append(i)
			db['names'][i] = name
	json.dump(db, fp, indent=2)
	
def query_database(query, db):
	qname, qset = shred(query)
	votes = dict()
	#find kmer matches in db
	for kmer in qset:
		name_ids = db['kmers'][kmer]
		for name_id in name_ids:
			if name_id not in votes:
				votes[name_id] = 0
			votes[name_id] = votes[name_id] + 1
			
	#figure out which db kmer set was matched to the most
	best_name_id =0 
	best_count = 0
	for name_id, count in votes.items():
		if count > best_count:
			best_name_id = name_id
			best_count = count
	actual_name = db['names'][str(best_name_id)]
	
	#return the name of the db reference sequence that the query kmers most matched to
	return actual_name
	



if __name__ == '__main__':
	print "building database"
	build_database(SeqIO.parse(open('four_reads.fasta', 'r'), 'fasta'), open('db.json', 'w'))
	print "querying database"
	print query_database(SeqIO.parse(open('one_read.fasta', 'r'), 'fasta').next(), json.load(open('db.json', 'r')))