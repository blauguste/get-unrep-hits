from Bio import SeqIO
from Bio import AlignIO
from Bio import Entrez
import pandas as pd
import sys
import pickle

Entrez.email = 'hdutcher@pdx.edu'

def get_species_name(nt_identifier):
    handle = Entrez.efetch(db='nucleotide', id=nt_identifier, retmode='xml')
    record = Entrez.read(handle)
    name = record[0]['GBSeq_organism']
    return name

def get_bottom_hits(bn, seed_alignment, restable, iter_ct):    

    seed_sto = AlignIO.read(open(seed_alignment, 'r'), 'stockholm')
    seed_accs = [a.id.split('/')[0] for a in seed_sto]
    
    if int(iter_ct) == 2:
        seed_names = [get_species_name(a) for a in seed_accs]
        seed_gs = [b.split(' ')[0] + ' ' + b.split(' ')[1] for b in seed_names]
    else:
        seed_gs = pickle.load(open(bn + '_seed_names.p', 'rb'))    

    cmsearch = ['target_name', 'target_accession', 'query_name', 'query_accession', \
            'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', \
            'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'genus', 'species']

    df = pd.read_csv(restable, delim_whitespace=True, header=None, names=cmsearch, \
                    skiprows=2, skipfooter=10, engine='python', index_col=False)

    df['gs'] = df['genus'] + ' ' + df['species']
    new_ss = pd.DataFrame()
    new_ss_names = []
    score = 0
    for i, row in df.sort_values(by=['E-value'], ascending=False).iterrows():
        model_name = row['query_accession']
        if row['gs'] in seed_gs:
            break
        elif row['gs'] in new_ss_names:
            continue
        elif score == row['score']:
            continue
        else:
            score = row['score']
            new_ss_names.append(row['gs'])
            new_ss = new_ss.append(row)

    if not new_ss.empty:
        seed_gs.extend(new_ss_names)
        new_ss['seq_from'] = new_ss['seq_from'].astype(int)
        new_ss['seq_to'] = new_ss['seq_to'].astype(int)

        pickle.dump(seed_gs, open(bn + '_seed_names.p', 'wb'))
        new_ss['seed_name'] = new_ss['target_name'] + '/' + new_ss['seq_from'].astype(str) + '-' + new_ss['seq_to'].astype(str)
        fn = bn + '_v' + str(iter_ct) + '.csv'
        new_ss[['seed_name', 'seq_from', 'seq_to', 'target_name']].to_csv(fn, sep='\t', header=None, index=None)

if __name__ == '__main__':
    if len(sys.argv) == 5:
        get_bottom_hits(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("Usage: get_bottom_hits.py basename seed_alignment cmsearch_results_table iteration_count")


