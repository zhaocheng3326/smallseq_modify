from __future__ import division
from collections import defaultdict
import sys 
#sys.path.insert(0, '/home/chenzh/PC/SM/small_RNA/core')
sys.path.insert(0, './')
import dr_tools, pysam
import argparse, os

class betweenRE(dr_tools.Cregion):
	pass


def bam_to_windows(inbam,out):
	inbamPysamObj = pysam.Samfile(inbam, "rb" )
	#p = inbam.split("/")
	#samplename = p[-2]
	tempCountfile = out + "_tmpCount.txt"
	finalCountfile = out + "_Count.txt"
	read2overlapCoords=defaultdict(list)

	for read in inbamPysamObj:
		readname = read.qname
		tid = read.rname
		readchr  = inbamPysamObj.getrname(tid)
		readstart = int(read.pos)
		readend = read.aend
		if read.is_reverse: 
			strand="-"
		else:
			strand="+"
		readlen = len(read.seq) #this is the actual read length (41M, means readlen=41)
		read_len = read.qlen  #this only considers matches (8S30M, means read_len=30)
 
		midpos = (readstart + readend)//2

		#retrieve list of overlapping coordinates
		overlap_list = betweenRE.overlappingpoint(readchr, midpos, strand)
		annotatedCount = len(overlap_list)

		#make a dictionary of read and overlapping coordinates
		read2overlapCoords[readname].append(overlap_list)

	with open(tempCountfile, "w") as outfh:
		for read in read2overlapCoords:
			coordsList = read2overlapCoords[read]
			readCount = len(coordsList)
			annotatedCount = readCount-coordsList.count([])
			#len(coordsList) is never zero
			for coord in coordsList:
				if len(coord) == 0:
					print(dr_tools.join(read, "NA", readCount, annotatedCount),file=outfh)
				else:
					###coord[1] will be double-counting
					coord = str(coord[0])  #otherwise I got keyError. it was "instance" type variable
					geneid = coord2geneid.get(coord, 'NA')
					print(dr_tools.join(read, geneid, readCount, annotatedCount),file=outfh)
	outfh.close()

	## readCount, annotatedCount scenarios
	# 1, 1  unique map, annotated to single gene, counts as 1
	# 2, 1  multi map, annotated to single gene, count as 1, discard other alignment
	# n, n  multi map, annotated to two genes, count as 1/n ?
	# k, m  where k>m and m>1, multi map, annotated to multi genes, count 1/m, discard other alignment
	#
	#formula is always: count = 1/annotatedCount
	geneid2counts={}
	unannotReadsDict={}
	for line in open(tempCountfile, "r"):
		p = line.split()
		read, geneid, readCount, annotatedCount = p
		annotatedCount = int(annotatedCount)
		if not geneid in geneid2counts: geneid2counts[geneid] = 0
		if annotatedCount > 0:
			geneid2counts[geneid] += 1/annotatedCount
		else:
			geneid2counts[geneid] += 0
		if annotatedCount < int(readCount) and annotatedCount == 0:
			unannotReadsDict[read] = 1

	num_unannot = len(unannotReadsDict)
	num_annot = 0
	for geneid in geneid2counts:
		if "P-cel" in geneid: continue
		if geneid == "NA": continue
		num_annot += geneid2counts[geneid]

	with open(finalCountfile, "w") as outfh2:
		print(dr_tools.join("#samples",os.path.basename(inbam).rstrip(".sort.bam")),file=outfh2)
		print(dr_tools.join("#unannotatedmolc", num_unannot),file=outfh2)
		print(dr_tools.join("#annotatedmolc", num_annot),file=outfh2)
		for geneid in geneidlist:
			print (dr_tools.join(geneid2name[geneid], geneid, geneid2counts.get(geneid, "0")),file=outfh2)
	outfh2.close()

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--inputbam', required=True)
	parser.add_argument('-o', '--output', required=True, help='directory for final counts')
	parser.add_argument('-g', '--GenePred', default="annotations/combined_annots.gp") 
	o = parser.parse_args()

	### read coordinates and add to windows
	geneid2name={}
	coord2geneid={}
	geneidlist=[]
	for line in open(o.GenePred, 'r'):
		p = line.split()
		chrom, startpos, endpos, strand, geneid, genename = p[0], int(p[1]), int(p[2]), p[5], p[3], p[6]
		#"add to windows
		r = betweenRE(chrom, startpos, endpos, strand)
		r.addtowindows()
		coord = chrom+":"+str(startpos+1)+"-"+str(endpos)+":"+strand
		coord2geneid[coord] = geneid #check length
		geneid2name[geneid] = genename
		geneidlist.append(geneid)

	### create outfolder, add sample names to list
	#for sample in sample_names:
		##prepare input files
	#	uniqmultibam = os.path.join(o.instardir, sample, "%s.bam" %sample)
	#	samplenames_with_fullpath.append(uniqmultibam)
	#	path_outbam = os.path.join(o.outdir, sample)

	### call function in parallel 
	bam_to_windows(o.inputbam,o.output)
