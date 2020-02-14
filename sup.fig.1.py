'''
For regenerating supplementary figure 1, 
the python code titled “sup.fig.1.py” is provided. 
The input files for this script are: 
1) ProxReg output files of all 90 Chip-Seq data, 
2) GO term related information extracted from org.Hs.eg.db. 
Once these two datasets were prepared, 
users can run the python script to generate the 
data for supplementary figure 1. 
The file paths in the python script need to be changes, 
users can follow the comments in the script to change related paths. 
'''

import sys,os

def get_go2gene():
	'''
	read the GO info
	'''
	out = {}
	filepath = "reverse_query_initial_GO_2.txt"
	with open(filepath) as f:
		for line in f:
			if line.startswith("Geneset.Type") is False and "Ontology Biological Process" in line:
				line = line.strip()
				goid = line.split("\t")[1]
				genelist= (line.split("\t")[3]).split(",")
				if len(genelist)<=1000:
					out[goid] = genelist
	return(out)


def returngotermfromfile(file):
	out = []
	with open(file) as f:
		for line in f:
			line = line.strip()
			if line.startswith("SYMBOL"):
				next
			else:
				go = line.split()[1]
				out.append(go)
	return(out)

def combinetwodict(tfgolist,gogenelist):
	inTF = {}
	notinTF = {}
	for k,v in tfgolist.items():
		inTF[k] = {}
		notinTF[k] = {}
		for kk,vv in gogenelist.items():
			if kk not in v:
				notinTF[k][kk] = gogenelist[kk]
	for k,v in out.items():
		print(k)
		for i,j in v.items():
			print(i+"\t"+str(len(j)))
	return(inTF,notinTF)


def get_gofortf():
	out = {}
	folder = ".../TF.from.Gold.Standard.v3/" # path to the GO data extracted from org.Hs.eg.db
	files  = os.listdir(folder)
	for i in files:
		name = i.split(".")[0]
		golist = returngotermfromfile(folder+i)
		out[name] = golist
	return(out)

def get_peaks_res_filelist():
	out = []
	path = ".../out/" # path ot the ProxReg results
	files = os.listdir(path)
	for i in files:
		if i.endswith("promoter_peaks.tab"):
			out.append(path+i)
	return(out)

def get_out_final(outlist,total):
	out = []
	tp = []
	for i in outlist:
		ids = i.split("\t")[1]
		tp.append(ids)
	name = (outlist[0].split("\t")[0]).split(".")[1]
	if name == "YY1SC281":
		name = "YY1"
	else:
		name = (outlist[0].split("\t")[0]).split(".")[1]
	cell = outlist[0].split("\t")[0]
	sub = total[name]
	for k,v in sub.items():
		get = 0
		for j in v:
			if j in tp:
				get = get+1
		out.append(cell+"\t"+k+"\t"+str(get)+"\t"+str(len(v)))
	return(out)


def remove_same_gene_from_go_group(dic):
	out={}
	for k,v in dic.items():
		tp = {}
		for kk,vv in v.items():
			for i in vv:
				tp[i] = 1
		out[k] = list(tp.keys())
	return(out)

def test_gene():
	dic={}
	file = ".../reverse_query_initial_GO_2.txt"
	with open(file) as f:
		for line in f:
			line = line.strip()
			if line.startswith("Geneset.Type") is False and "Ontology Biological Process" in line:
				geneset = line.split("\t")[3]
				tp = geneset.split(",")
				for i in tp:
					dic[i] = 1
	return(list(dic.keys()))



def get_final(files,dic,fck):
	print("write results into file...")
	path = ".../out.data/" #path to the coutn files
	outfile=open(".../test.txt",'w') #output file path
	outfile.write("data\tgene.peak.targetGO\ttotalgenetargetGO\tgene.peak.notargetgene\ttotalgenenotargetgene\n")
	for i in files:
		if "Pol2" in i:
			next
		else:
			tp = {}
			sp = {}
			base = os.path.basename(i)
			name = ".".join(base.split(".")[0:2])
			name = name.upper()
			if "YY1SC281" in name:
				key = "YY1"
			else:
				key  = name.split(".")[1]
			with open(i) as f:
				for line in f:
					line = line.strip()
					if line.startswith("peak_id"):
						next
					else:
						gene_id = line.split("\t")[4]
						if "e" in line.split("\t")[9]:
							dtss = 0
						else:
							dtss    = int(line.split("\t")[9])
						if abs(dtss)<=2000 and gene_id in dic[key]:
						 	tp[gene_id] = 1
						elif abs(dtss)<=2000 and gene_id in fck[key]:
							sp[gene_id] = 1
		outfile.write(name+"\t"+str(len(tp.keys()))+"\t"+str(len(dic[key]))+"\t"+str(len(sp.keys()))+"\t"+str(len(fck[key]))+"\n")
			
def get_not_in_genelist(intfgene,genelist):
	
	out = {}
	for k,v in intfgene.items():
		nots = list(set(v)^set(genelist))
		out[k] = nots
	return(out)

def get_promoter_geneid_countforonepeakfile(files,inTF):
	print("write results into file...")
	path = ".../out.data/" #path to a folder to save the count files
	for i in files:
		if "Pol2" in i:
			next
		else:
			tp = []
			base = os.path.basename(i)
			name = ".".join(base.split(".")[0:2])
			name = name.upper()
			with open(i) as f:
				for line in f:
					line = line.strip()
					if line.startswith("peak_id"):
						next
					else:
						gene_id = line.split("\t")[4]
						if "e" in line.split("\t")[9]:
							dtss = 0
						else:
							dtss    = int(line.split("\t")[9])
						if abs(dtss)<=2000:
						 	tp.append(name+"\t"+gene_id+"\t"+str(dtss))
			outlist = get_out_final(tp,inTF)
			outfile = open(path+name+".inTF.count.txt",'w')
			outfile.write("cell.TF\tGO\tNum.of.gene.peak.promoter\tNum.of.gene.GO.term\n")
			for j in outlist:
				outfile.write(j+"\n")
	print("Done!")

def get_promoter_geneid_countforonepeakfile_notinTF(files,notinTF):
	print("write results into file...")
	path = ".../out.data/"#path to a folder to save the count files
	for i in files:
		if "Pol2" in i:
			next
		else:
			tp = []
			base = os.path.basename(i)
			name = ".".join(base.split(".")[0:2])
			name = name.upper()
			with open(i) as f:
				for line in f:
					line = line.strip()
					if line.startswith("peak_id"):
						next
					else:
						gene_id = line.split("\t")[4]
						if "e" in line.split("\t")[9]:
							dtss = 0
						else:
							dtss    = int(line.split("\t")[9])
						if abs(dtss)<=2000:
						 	tp.append(name+"\t"+gene_id+"\t"+str(dtss))
			outlist = get_out_final(tp,notinTF)
			outfile = open(path+name+".notinTF.count.txt",'w')
			outfile.write("cell.TF\tGO\tNum.of.gene.peak.promoter\tNum.of.gene.GO.term\n")
			for j in outlist:
				outfile.write(j+"\n")
	print("Done!")

def main():
	gogenelist=get_go2gene()
	tfgolist=get_gofortf()
	files = get_peaks_res_filelist()
	genelist=test_gene()
	inTF,notinTF = combinetwodict(tfgolist,gogenelist)
	get_promoter_geneid_countforonepeakfile(files,inTF)
	get_promoter_geneid_countforonepeakfile_notinTF(files,notinTF)
	intfgene=remove_same_gene_from_go_group(inTF)
	notintfgene=get_not_in_genelist(intfgene,genelist)
	get_final(files,intfgene,notintfgene)

if __name__ == "__main__":
	main()









