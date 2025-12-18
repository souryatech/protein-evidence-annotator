import mysql.connector

geneid_details = {}

htab_file = "hmmscan.htab"
fasta_file = "prodigal2fasta.nostars.faa"
tmhmm_file = "prodigal2fasta.nostars.tmhmm.short"


# Parsing through fasta table

with open(fasta_file, 'r') as fasta:
    for line in fasta:
        if "1_" in line:
            geneid_details.update({line.replace("_polypeptide\n", "").replace(">", ""):{'hmm':None, 'blast':None, 'tmhmm':None}})



conn = mysql.connector.connect(
    user='root',
    host='localhost',
    password = 'Likeblue#55',
    database='annot'
) 

# Parsing through hmmscan table

with open(htab_file, 'r') as htab:
    for line in htab:
        row = line.split('\t')
        evalue = float(row[19])
        if evalue < 1e-50:
            key = row[5].replace("_polypeptide", "")
            if not geneid_details[key]['hmm']:
                geneid_details[key]['hmm'] = (row[15],evalue)
       
# import itertools
# print(dict(itertools.islice(geneid_details.items(), 0, 4)))
cursor = conn.cursor()

# Parsing through BLAST table

cursor.execute("""
    SELECT qry_id, product, evalue
    FROM blast
    WHERE evalue < 1e-50
""")
for qry_id, product, evalue in cursor.fetchall():
    if not geneid_details[qry_id]['blast']:
        geneid_details[qry_id]['blast'] = (product,evalue)
conn.commit()
cursor.close()

# print(dict(itertools.islice(geneid_details.items(), 0, 10)))


# Parsing through tmhmm scan table
with open(tmhmm_file, 'r') as tmhmm:
    for line in tmhmm:
        row = line.split('\t')
        if int(row[2][8:]) > 0:
            key = row[0].replace("_polypeptide", "")
            geneid_details[key]['tmhmm'] = True

# print(dict(itertools.islice(geneid_details.items(), 0, 10)))

#Output loop & file creation

with open('protein_evidence.txt', 'w') as output:
    for gene_id in geneid_details:
        if geneid_details[gene_id]['hmm']:
            str = f"{gene_id}\t{geneid_details[gene_id]['hmm'][0]}"
        elif geneid_details[gene_id]['blast']:
            str = f"{gene_id}\t{geneid_details[gene_id]['blast'][0]}"
        elif geneid_details[gene_id]['tmhmm']:
            str = f"{gene_id}\tPutative transmembrane protein"
        else:
            str = f"{gene_id}\tHypothetical Protein"
        print(str)
        output.write(str+"\n")
