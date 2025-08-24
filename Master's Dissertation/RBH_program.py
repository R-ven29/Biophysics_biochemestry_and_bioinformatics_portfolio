import re 
import Bio
import os
import pandas as pd
import subprocess
import dask.dataframe as dd
import time 
import concurrent.futures
from tqdm import tqdm
from unipressed import IdMappingClient
from unipressed import UniprotkbClient
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict

start_time = time.perf_counter()

#QUERY FASTA RETRIEVAL
Entrez.email="example@example.es"
IDs = input('protein ncbi ids separated by comas:')
with Entrez.efetch(db="protein", id=IDs, rettype="gb", retmode="text") as handle:
    SeqIO.write((record for record in SeqIO.parse(handle, "gb") if len(record.seq) >= 1),"queries.fasta","fasta")
    
#1ST SEARCH
if not os.path.exists('uniprot_flavoproteins_with_function.dmnd'):
    subprocess.run('wsl diamond makedb --in "uniprot_flavoproteins_with_function.fasta" --db "uniprot_flavoproteins_with_function"\
                   --no-parse-seqids')
            
subprocess.run('wsl diamond blastp --db "uniprot_flavoproteins_with_function.dmnd" --query queries.fasta\
               --outfmt 6 qseqid sseqid qlen slen pident evalue bitscore score length -k 10000 -o "out.txt" --ultra-sensitive')
    
out= pd.read_table("out.txt", names=['qseqid', 'sseqid', 'qlen', 'slen', 'pident', 'evalue', 'bitscore', 'score', 'length'])
 
#DATABASE FASTA DESCRIPTION SAVING 
Flavoprotein_ids=[]
records=SeqIO.parse("uniprot_flavoproteins_with_function.fasta", "fasta")
for record in records:
    Flavoprotein_ids.append(record.description)

#TAXONOMY
seqid_to_sciname = {}
for item in Flavoprotein_ids:
    sseqid_match = re.search(r"^(\S+)", item)  
    if sseqid_match:
        sseqid = sseqid_match.group(1)
        description_match = re.search(r"OS=(.+?)(?: OX=|\()", item)
        if description_match:
            seqid_to_sciname[sseqid] = description_match.group(1)
out['sciname'] = out['sseqid'].map(seqid_to_sciname)

#COVERAGE CRITERIA AND BEST HITS FILTERING
ddf = dd.from_pandas(out, npartitions=8)
out = ddf[((ddf['length'] / ddf['qlen'])> 0.4) & ((ddf['length'] / ddf['slen']) > 0.4)].compute()
best_indices = out.groupby('sciname')['evalue'].idxmin()
Best_hits_df = out.loc[best_indices].reset_index(drop=True)

#FROM NCBI QUERY CODE TO UNIPROT's
request = IdMappingClient.submit(source="EMBL-GenBank-DDBJ_CDS", dest="UniProtKB", ids=IDs.split(','))
results = []
for result in request.each_result():
    results.append(result)
    time.sleep(0.05)   
mapid=dict()
for i in range(len(results)):
    mapid[results[i]['from']]=results[i]['to']
Best_hits_df['qseqid_uniprot'] = Best_hits_df['qseqid'].map(mapid)

#ADDING PROTEIN DESCRIPTION TO UBH FROM FASTA 
seqid_to_description = {}
for item in Flavoprotein_ids:
    sseqid_match = re.search(r"^(\S+)", item)  
    if sseqid_match:
        sseqid = sseqid_match.group(1)
        description_match = re.search(r"(?<=\S) (.*?) (?=OS=)", item)
        if description_match:
            seqid_to_description[sseqid] = description_match.group(0)
Best_hits_df['description'] = Best_hits_df['sseqid'].map(seqid_to_description)

#RENAMING COLUMNS FOR BETTER OUTPUT
Best_hits_df.rename(columns= { 'qseqid_uniprot': 'query id UniProt', 
                                'qseqid': 'query id NCBI',
                                'sseqid' : 'subject id',
                                'qlen': 'query length',
                                'slen': 'subject length',
                                'pident': 'indentity percentage',
                                'evalue' : 'E-value',
                                'bitscore': 'bit score',
                                'score': 'score',
                                'length': 'alignment length',
                                'sciname': 'scientific name'}, inplace=True)

#ADDING PROTEIN INFORMATION TO UBH FROM UNIPROT
def fetch_uniprot_record(uniprot_id):
    try:
        data = UniprotkbClient.fetch_one(uniprot_id)
        result = {'EC 1': None, 'EC 2': None, 'pathway 1': None, 'pathway 2': None,'estechiometry 1': None,
                  'estechiometry 2': None, 'estechiometry 3': None, 'estechiometry 4': None, 'cofactor 1': None, 'cofactor 2': None,
                  'cofactor 3': None, 'cofactor 4': None, 'activity 1': None,'activity 2': None, 'RHEA 1': None, 'RHEA 2': None}
        pathways=[]
        activities=[]
        reactions=[]
        ec_numbers=[]
        cofactors=[]
        estechiometry=[]
        if 'comments' in data:
            comments = data['comments']     
            for comment in comments:
                if comment['commentType'] == 'COFACTOR':
                    if 'note' in comment:
                        if 'texts' in comment['note']:
                            estechiometry.append(comment['note']['texts'][0].get('value', None))
                            if len(estechiometry) > 1:
                                result['estechiometry 1'] = estechiometry[0]
                                result['estechiometry 2'] = estechiometry[1]
                            if len(estechiometry) > 2:
                                result['estechiometry 3'] = estechiometry[2]    
                            if len(estechiometry) > 3:
                                result['estechiometry 4'] = estechiometry[3]
                            if len(estechiometry) == 1:
                                result['estechiometry 1'] = estechiometry[0]
                                result['estechiometry 2'] = None
                                result['estechiometry 3'] = None
                                result['estechiometry 4'] = None
                    if 'cofactors' in comment:
                        cofactors.append(comment['cofactors'][0]['name'])
                        if len(cofactors) > 1:
                            result['cofactor 1'] = cofactors[0]
                            result['cofactor 2'] = cofactors[1]
                        if len(cofactors) > 2:
                            result['cofactor 3'] = cofactors[2]           
                        if len(cofactors) > 3:
                            result['cofactor 4'] = cofactors[3]
                        if len(cofactors) == 1:
                            result['cofactor 1'] = cofactors[0]
                            result['cofactor 2'] = None
                            result['cofactor 3'] = None
                            result['cofactor 4'] = None       
                if comment['commentType'] == 'PATHWAY':
                    pathways.append(comment['texts'][0].get('value', None))
                if len(pathways) > 1:
                    result['pathway 1'] = pathways[0]
                    result['pathway 2'] = pathways[1]
                if len(pathways) == 1:
                    result['pathway 1'] = pathways[0]
                    result['pathway 2'] = None
                if comment['commentType'] == 'CATALYTIC ACTIVITY' and 'reaction' in comment:
                    reaction =comment['reaction'] 
                    if 'name' in reaction :
                        activities.append(reaction.get('name', None)) 
                    if len(activities) > 1:
                        result['activity 1'] = activities[0]
                        result['activity 2'] = activities[1]
                    if len(activities) == 1:
                        result['activity 1'] = activities[0]
                        result['activity 2'] = None                                               
                    if 'reactionCrossReferences' in reaction and len(reaction['reactionCrossReferences']) > 0:
                        reactions.append(reaction['reactionCrossReferences'][0].get('id', None))
                    if len(reactions) > 1:
                        result['RHEA 1'] = reactions[0]
                        result['RHEA 2'] = reactions[1]
                    if len(reactions) == 1:
                        result['RHEA 1'] = reactions[0]
                        result['RHEA 2'] = None
                    if 'ecNumber' in reaction:
                        ec_numbers.append(reaction['ecNumber'])
                    if len(ec_numbers) > 1:
                        result['EC 1'] = ec_numbers[0]
                        result['EC 2'] = ec_numbers[1] 
                    if len(ec_numbers) == 1:
                        result['EC 1'] = ec_numbers[0]
                        result['EC 2'] = None     
        return result
    except Exception as e:
        return {'EC 1': None, 'EC 2': None, 'pathway 1': None, 'pathway 2': None,'estechiometry 1': None,
                'estechiometry 2': None, 'estechiometry 3': None, 'estechiometry 4': None, 'cofactor 1': None, 'cofactor 2': None,
                'cofactor 3': None, 'cofactor 4': None, 'activity 1': None,'activity 2': None, 'RHEA 1': None, 'RHEA 2': None}

Best_hits_df['sseqid_short'] = Best_hits_df['subject id'].str.extract(r"\|([^|]+)\|").astype(str)
sseqid_uniprot=list(Best_hits_df['sseqid_short'])

with concurrent.futures.ThreadPoolExecutor() as executor:
    results_UBH = list(tqdm(executor.map(fetch_uniprot_record, sseqid_uniprot),
            total=len(sseqid_uniprot),
            desc="Fetching UniProt records"))

UBH_annotations=pd.DataFrame(results_UBH)
Best_hits_df = pd.concat([Best_hits_df, UBH_annotations], axis=1)

#SAVING UBH RESULTS
Best_hits_df.to_excel('UBH_out.xlsx', index=False)

#MAKING 2ND SEARCH QUERY FASTA FILE WITH UBH
with open(r'uniprot_flavoproteins_with_function.fasta', 'r') as myco_r:
    sequences = defaultdict(list)
    b =''
    for line in myco_r:
        line = line.strip()
        if line.startswith('>'):
            sequences['header'].append(line)
            if b !='':
                sequences['sequence'].append(b)
            b=''
        else:
            b += line.strip()
    if b!='':
            sequences['sequence'].append(b)
            
sequences['header'] = [re.search(r'>(\S+)', header).group(1) for header in sequences['header']]
mapseq = {header: sequence for header, sequence in zip(sequences['header'], sequences['sequence'])}

BH_sequences=Best_hits_df['subject id'].map(mapseq)
with open('best_hits.fasta', 'w') as myco_w:
    for sseqid, sequence in zip(Best_hits_df['subject id'], BH_sequences):
        if pd.notna(sequence): 
            myco_w.write(f">{sseqid}\n")
            myco_w.write(f"{sequence}\n")

#2ND SEARCH
if not os.path.exists('uniprotkb_Mycobacterium_tuberculosis_H37Rv.dmnd'):
    subprocess.run('wsl diamond makedb --in "uniprotkb_Mycobacterium_tuberculosis_H37Rv.fasta"\
                   --db "uniprotkb_Mycobacterium_tuberculosis_H37Rv" --no-parse-seqids')
            
subprocess.run('wsl diamond blastp --db "uniprotkb_Mycobacterium_tuberculosis_H37Rv.dmnd" --query best_hits.fasta\
               --outfmt 6 qseqid sseqid qlen slen pident evalue bitscore score length -k 1000 -o "out_reciprocal.txt" --ultra-sensitive')

out_reciprocal= pd.read_table("out_reciprocal.txt", names=['qseqid', 'sseqid', 'qlen', 'slen', 'pident', 'evalue', 'bitscore', 'score', 'length'])    

#COVERAGE CRITERIA AND BEST HITS FILTERING
ddf_r = dd.from_pandas(out_reciprocal, npartitions=8)
out_reciprocal = ddf_r[((ddf_r['length'] / ddf_r['qlen']) > 0.4) & ((ddf_r['length'] / ddf_r['slen']) > 0.4)].compute()
out_reciprocal['sciname'] = out_reciprocal['qseqid'].map(seqid_to_sciname)
best_indices_r = out_reciprocal.groupby('sciname')['evalue'].idxmin()

#MERGING TABLES TO DISCOVER RECIPROCAL HITS
out_reciprocal = out_reciprocal.loc[best_indices_r].reset_index(drop=True)

out_reciprocal['sseqid_uniprot'] = out_reciprocal['sseqid'].str.extract(r"\|([^|]+)\|")
out_reciprocal['sseqid_uniprot'] = out_reciprocal['sseqid_uniprot'].astype(str)

Best_hits_df['query id UniProt'] = Best_hits_df['query id UniProt'].astype(str)

best_hits_out = pd.merge(out_reciprocal, Best_hits_df, how='inner', left_on=['sseqid_uniprot', 'qseqid'], right_on=['query id UniProt', 'subject id'])
best_hits_out.drop(columns=['sseqid_uniprot', 'qseqid', 'sseqid', 'qlen', 'slen', 'pident', 'evalue', 'bitscore', 'score_y', 'length'], inplace=True)
best_hits_out.drop_duplicates(inplace=True)

#ADDING DESCRIPTION FOR RBH FROM FASTA
best_hits_out['description'] = best_hits_out['subject id'].map(seqid_to_description)

#RENAMING COLUMNS FOR BETTER OUTPUT
best_hits_out.rename(columns= {'score_x': 'score'}, inplace=True)
best_hits_out['subject id'] = best_hits_out['subject id'].str.extract(r"\|([^|]+)\|").astype(str)

#SAVING RBH RESULTS
best_hits_out.to_excel('RBH_out.xlsx', index=False)

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.6f} seconds")