import pandas as pd
import tempfile
import subprocess
from io import StringIO
import re
from collections import defaultdict
import os
import argparse
import datetime
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-p', '--profile',
        type=str,
        help='Input profile(s)',
        required=True,
        dest='profile',
        nargs='+'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        help='Output folder to store CAMI formats in',
        required=True,
        dest='output'
    )


    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=1,
        help='number of threads to do multiprocessing on',
        dest='threads',
    )
    parser.add_argument(
        '-c', '--check',
        action='store_true',
        help='Check if output file already exists',
        dest='check'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='be verbose',
        dest='verbose'
    )

    return parser.parse_args()

    
p_qids = {
    'Silva_20200116': re.compile(r'^(\S+)[\s\;\.]') ,
    'genomic_20220524': re.compile(r'^(\S+)[\s\;\.]') ,
    'db_mOTU_20221205': re.compile(r'^(\S+)[\s\;\.]') ,
    'metaphlan_20221125': re.compile(r'^(\S+)[\s\;\.]') ,
    'bac120_marker_genes_r207_k16_20220726': re.compile(r'^(\S+\.\d+)\_')  
}

#ranks = ['superkingdom','phylum','class','order','family','genus','species']

def get_sampleinfo(profile, format):

    p_kma = re.compile(r'(\w+)\_(genomic|Silva|mOTUs|metaphlan|gtdb)\.mapstat')
    p_motus = re.compile(r'taxonomy\_(\w+)\_(profile|BIOM)')
    p_bracken = re.compile(r'(\w+)_(species|genus|family|order|class|combined)\.bracken\.')

    if format == 'kma':
        m = p_kma.findall(profile)[0]
        sampleid = m[0]

        output = "_".join(m) + '_KMA.cami'
    elif format == 'motus':
        m = p_motus.findall(profile)[0]
        sampleid = m[0]
        output = sampleid + "_mOTUs.cami"
    elif format == "bracken":
        m = p_bracken.findall(profile)[0]
        sampleid = m[0]
        output = "_".join(m) + '_bracken.cami'


    return sampleid, output

def generate_header(sample, ranks):
    CAMI_header = "@SampleID:{sample}\n@Version:0.9.3\n@Ranks:{ranks}\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n"
    return CAMI_header.format(sample=sample, ranks=ranks)


def get_taxa(df, database):
    if database.startswith('metaphlan'):
        query = "use taxonomy; select * from {} where {};".format(
            database,
            " OR ".join(df['qid'].dropna().apply(lambda x: f"id like '{x}%'", meta=('qid', 'str'))
        )
    else:
        query = "use taxonomy; select * from {} where superkingdom_tax is not null and id in ({});".format(
            database,
            ",".join(df['qid'].compute().dropna().apply(lambda x: "\'" + x + "\'").tolist())
        )    
    with tempfile.NamedTemporaryFile(mode='w') as f:
        f.write(query)
        f.flush()

        p = subprocess.run(f"mysql < {f.name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if p.returncode == 0:
            taxInfo = pd.read_csv(StringIO(p.stdout.decode()), sep='\t', low_memory=False)
            taxInfo.drop_duplicates(inplace=True)
            taxInfo['id_len'].fillna(taxInfo['genome_size'], inplace=True)
            taxInfo.dropna(subset=['id_len'], inplace=True)
            #if taxInfo['id_len'].isna().all():
            #    taxInfo.drop(columns='id_len', inplace=True)
            #    mref_gene_lengths = ref_gene_lengths.loc[ref_gene_lengths['db'] == database, ].copy()
            #    mref_gene_lengths['qid'] = mref_gene_lengths['id'].str.extract(p_qids[database])
            #    if taxInfo['id'].isin(mref_gene_lengths['id']).any():
            #        taxInfo = taxInfo.merge(mref_gene_lengths[['id', 'length']], on='id')
            #    elif taxInfo['id'].isin(mref_gene_lengths['#name']).any():
            #        taxInfo = taxInfo.merge(mref_gene_lengths[['#name', 'length']], left_on='id', right_on='#name')
            #    elif taxInfo['id'].isin(mref_gene_lengths['qid']).any():
            #        taxInfo = taxInfo.merge(mref_gene_lengths[['qid', 'length']], left_on='id', right_on='qid')
            #    else:
            #        taxInfo['id'] = taxInfo['id'].str.extract(p_qids[database])
            #        taxInfo = taxInfo.merge(mref_gene_lengths[['id', 'length']], on='id')

             #   taxInfo.rename(columns={'length': 'id_len'}, inplace=True)
            return taxInfo
        else:
            print(p.stderr.decode(), file=log)



def build_paths(ranks, taxInfo):
    
    # taxid, rank, taxpath, taxpathsn
    camiTaxa = []
    prev_ranks = []
    for rank in ranks:
        matching_cols = taxInfo.columns[taxInfo.columns.str.contains(rank +'_' )].tolist()
        prev_ranks += matching_cols

        matching_rows = taxInfo[matching_cols + ['id', 'id_len']].drop_duplicates()#.dropna(subset=matching_cols)
        for rowid, row in matching_rows.iterrows():
            if taxInfo.loc[rowid, ['id'] + prev_ranks[0::2]].isna().any():
                continue
                
            taxid = int(row[matching_cols[1]])
            qid = row['id']
            
            taxpathsn = "|".join(taxInfo.loc[rowid, prev_ranks[0::2]].fillna('').values.tolist())
            taxpath = "|".join(taxInfo.loc[rowid, prev_ranks[1::2]].fillna('').astype(str).values.tolist())
            
            d = {}
            d['TAXID'] = taxid
            d['RANK'] = rank
            d['TAXPATHSN'] = taxpathsn
            d['TAXPATH']  = taxpath
            d = pd.DataFrame.from_dict(d, orient='index').T
            d['qid'] = qid
            d['qid_len'] = row['id_len']
            
            camiTaxa.append(d)

    return pd.concat(camiTaxa, axis=0)

def gtdb_get_id(x):
    p = re.compile(r'^\w{2}\_(\S+\.\d+)\_') 
    try:
        m = p.search(x)
        return m.group(1)
    except AttributeError:
        return pd.NA

def get_database(mapstatfile):
    # ! grep '## database' 'data/plant_associated_sample9_Silva.mapstat' | cut -f2
    cmd = f"grep '## database' {mapstatfile} | cut -f2"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return p.stdout.decode().strip()

taxdumpDir = '/home/databases/metagenomics/db/gtdb_bac120_20220726/gtdb-taxdump/R207/'
gtdbdumpDir = '/home/databases/metagenomics/db/gtdb_bac120_20220726/'

meta_lineages = pd.DataFrame(
    {x: 'object' for x in ['superkingdom_name', 'superkingdom_tax', 'phylum_name', 'phylum_tax',  'class_name', 'class_tax', 'order_name', 'order_tax', 'family_name','family_tax', 'genus_name', 'genus_tax', 'species_name','species_tax']},
    index=[0]
)

def gtdb_get_lengths(ref):
    p = subprocess.run(f"grep {ref} {os.path.join(gtdbdumpDir, 'bac120_marker_genes_reps_r207.seqkit.length')}", shell=True, stdout=subprocess.PIPE)
    if p.returncode == 0:
        o = p.stdout.decode().strip().split()
        return int(o[-1])
    else:
        print(p)
        return pd.NA

def gtdb_get_taxa(df):
    taxInfo = df[['qid', '# refSequence']].copy()
    taxInfo['id'] = taxInfo['# refSequence'].apply(gtdb_get_id, meta=('id', 'str'))#.compute()
    taxInfo['taxid'] = taxInfo['id'].apply(gtdb_get_taxid, meta=('taxid', 'int'))
    taxInfo['id_len'] = taxInfo['# refSequence'].apply(gtdb_get_lengths, meta=('id_len', 'str')) 
    
    lineages = taxInfo.apply(lambda x: gtdb_get_lineage(x['taxid']), result_type='expand', axis=1, meta=meta_lineages)
    
    full = dd.concat([taxInfo, lineages], axis=1, interleave_partitions=True)
    #full.rename(columns={'#name': 'id', 'length': 'id_len', 'taxid': 'strain_tax'})
    db='bac120_marker_genes_r207_k16_20220726'
    full['kma_db_name'] = db
    full = full.drop(columns=['id'])
    full = full.rename(columns={'taxid': 'strain_tax', 'qid': 'id'})
    with ProgressBar():
        full = full.compute()
    return full

def gtdb_get_taxid(id):
    
    cmd = f"grep {id} {taxdumpDir}/taxid.map"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode == 0:
        o = p.stdout.decode().strip()
        taxid = o.split('\t')[-1]
    else: 
        taxid = pd.NA
    return taxid


def gtdb_get_lineage(taxid):
    cmd = f"echo {taxid} | taxonkit lineage --data-dir {taxdumpDir} -t -R"
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o = p.stdout.decode().strip()
    lineages = o.split('\t')[1:]
    lineages = [l.split(';') for l in lineages]
    lineages = gtdb_convert_to_db_format(lineages)
    return lineages #lineage_names, lineage_ids, lineage_ranks

def gtdb_convert_to_db_format(lineages):
    result = pd.Series()
    for taxname, taxid, taxrank in zip(*lineages):
        if taxrank == 'no rank': continue
        result[taxrank + '_name'] = taxname 
        result[taxrank + '_tax'] = taxid 
    return result

def mapstat2cami(mapstatfile, outfile, ranks, verbose=False):

    # figure out which reference database was used
    database = get_database(mapstatfile)
    
    # load mapstat file
    df = pd.read_csv(mapstatfile, sep='\t', skiprows=6)
    df = df.loc[df['fragmentCountAln'] > 0]
    df['qid'] = df['# refSequence'].str.extract(p_qids[database])
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: loaded profile ({df.shape})", file=log)
    ddf = dd.from_pandas(df, npartitions=20)
    #ddf = ddf.repartition(partition_size='100MB')
    
    # get taxonomic information for ref sequences
    
    if database.startswith('bac120'):
       taxInfo = gtdb_get_taxa(ddf) 
    else:
        taxInfo = get_taxa(ddf, database)
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: loaded taxonomic info for {database} ({taxInfo.shape})", file=log)
    
    # build taxonomic paths according to cami format
    taxaPaths = build_paths(ranks, taxInfo)
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: Built CAMI tax paths", file=log)
    
    
    # merge
    dfMerged = taxaPaths.merge(df[['qid', 'fragmentCountAln']], on='qid', how='left')
    if dfMerged['fragmentCountAln'].isna().all():
        taxaPaths['qid'] = taxaPaths['qid'].str.extract(p_qids[database])
        dfMerged = taxaPaths.merge(df[['qid', 'fragmentCountAln']], on='qid', how='left')
    # local clade abundances
    dfMerged['local_abundance'] = dfMerged['fragmentCountAln'] / (dfMerged['qid_len'] / 1e3)
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: Calculated local clade abundances", file=log)
    
    #  sum by clade and calculate clade_read_abundances (not relative)
    dfClades = dfMerged.groupby(['TAXID', 'RANK', 'TAXPATHSN', 'TAXPATH']).agg({'local_abundance': 'sum', 'qid_len': 'mean'})
    dfClades['clade_abn'] = dfClades['local_abundance'] * (dfClades['qid_len']/1e3)
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: Calculated clade read abundances", file=log)
    
    # calculcate sum clade abundances
    sumCladeAbundances = dfClades.reset_index().groupby(['RANK']).agg({'clade_abn': 'sum'})
    sumCladeAbundances.rename(columns={'clade_abn': 'sum_clade'}, inplace=True)
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: Calculated total clade abundances", file=log)
    
    # merge into clade sums
    dfClades2 = dfClades.reset_index().merge(sumCladeAbundances, on='RANK', how='left')
    
    # calculate clade relative abundances
    dfClades2['PERCENTAGE'] = (dfClades2['local_abundance'] / dfClades2['sum_clade']) * 100
    if verbose: 
        with open(logfile, 'a') as log:
            print(f"{mapstatfile}: Calculcated relative abundances", file=log)
    
    CAMIdata = dfClades2[['TAXID', 'RANK', 'TAXPATH', 'TAXPATHSN', 'PERCENTAGE']].copy()
    CAMIdata = CAMIdata.groupby(['TAXID', 'RANK', 'TAXPATH', 'TAXPATHSN']).agg({'PERCENTAGE': 'sum'})
    CAMIdata.rename(columns={'TAXID':'@@TAXID' }, inplace=True)
    if CAMIdata['PERCENTAGE'].isna().all():
        if verbose:
            with open(logfile, 'a') as log:
                print(f"Failed to write to CAMI output file: {outfile}.. cleaning up", file=log)
        subprocess.run(f"rm {outfile}", shell=True)
    else:
        CAMIdata.sort_values(
            by="RANK", 
            key=lambda column: column.map(lambda e: ranks.index(e))
        ).to_csv(outfile, sep='\t', mode='a', header=None)
        if verbose: 
            with open(logfile, 'a') as log:
                print(f"{mapstatfile}: Written to CAMI output file: {outfile}", file=log)

def kma2cami(profile, ranks, outdir='.', verbose=False):
    sampleid, output = get_sampleinfo(profile, 'kma')
    
    outfile = os.path.join(outdir, output)
    header = generate_header(sample=sampleid, ranks="|".join(ranks))
    with open(outfile, 'w') as outstream:
        outstream.write(header)
    
    mapstat2cami(mapstatfile=profile, outfile = outfile, ranks=ranks, verbose=verbose)
    
  

if __name__ == "__main__":
    args = parse_args()

    # create output
    os.makedirs(args.output, exist_ok=True)

    global logfile
    logfile = os.path.join(args.output, 'kma2cami.log')
    if args.verbose:
        log = open(logfile, 'a')
        print(f"## {datetime.datetime.now()} ##", file=log)
        log.close()

    # define ranks
    ranks = ['superkingdom','phylum','class','order','family','genus','species']
    
    #global ref_gene_lengths
    #ref_gene_lengths = pd.read_csv('all_gene_lengths.csv', low_memory=False)

    for profile in args.profile:
        kma2cami(profile, outdir=args.output, verbose=args.verbose, ranks=ranks)
