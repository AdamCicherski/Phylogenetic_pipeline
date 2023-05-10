import requests
from Bio import SeqIO
from io import StringIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Phylo
import subprocess
import time
import os
import re
import logging
import argparse
import shutil
from requests.adapters import HTTPAdapter, Retry

def exec_path(p):
    if shutil.which(p) is not None:
        return p
    else:
        raise FileNotFoundError(p)

def console_args():
#Parsing command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=exec_path, default='iqtree2', help='path to iqtree')
    parser.add_argument('-cd', type=exec_path, default='cd-hit', help='path to cd-hit')
    parser.add_argument('-s', type=exec_path, default='./clann', help='path to clann')
    parser.add_argument('--cl', type=exec_path, default='clustalo', help='path to clustalo')
    parser.add_argument('-n', type=int, default=1, help='number of threads')
    parser.add_argument('-o', type=str, default='phylo_results', help='path to output directory')
    parser.add_argument('--og', type=str, default='Plicaturopsis_crispa', help='outgroup name')
    parser.add_argument('-t', type=str, default='boletales', help='taxon name')
    parser.add_argument('-D', action='store_true')


    return parser

# Setting up logger
def set_logger():
    logger = logging.getLogger('Phylo_aplication')
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    fileHandler = logging.FileHandler("phylo_pipeline.log", 'w')
    fileHandler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fileHandler.setFormatter(formatter)
    formatter2 = logging.Formatter('%(levelname)s - %(message)s')
    console_handler.setFormatter(formatter2)
    logger.addHandler(console_handler)
    logger.addHandler(fileHandler)
    return logger








re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        yield response
        batch_url = get_next_link(response.headers)


def get_fasta(url):
    s=[]
    for batch in get_batch(url):
        s += list(SeqIO.parse(StringIO(batch.text), 'fasta'))
    return s



def get_organisms_names():
    pass

def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        kw.get('logger', logging.getLogger(__name__)).info('function {} performance {} sec '.format(method.__name__, (te - ts).__round__(2)))
        return result

    return timed

# Download proteomes sequences from uniprot database
@timeit
def download_data(taxa, outgroup, paths, logger):

    url = 'https://rest.uniprot.org/proteomes/stream?&format=list&query=%28{}%29'.format(outgroup)
    try:
        outgroup_res = requests.get(url).text.split('\n')[0]
        if not re.match('UP\d+', outgroup_res):
            raise ValueError('No proteome for given out group name in Uniprot proteomes data base')
    except requests.exceptions.RequestException as e:
        logger.error(e)
        raise SystemExit(e)
    except ValueError as e:
        logger.error(e)
        SystemExit(e)

    url = 'https://rest.uniprot.org/proteomes/stream?&format=list&query=%28{}%29'.format(taxa)

    try:
        r = requests.get(url)
        results = []
        proteomes = r.text.split() + [outgroup_res]
        print(proteomes)
    except requests.exceptions.RequestException as e:
        logger.error(e)
        raise SystemExit(e)

    for i, proteome in enumerate(proteomes):
        print(proteome)
        url2 = 'https://rest.uniprot.org/uniprotkb/search?&format=fasta&query=%28%28proteome%3A{}%29%29&size=500'.format(proteome)
        try:
            s =  get_fasta(url2)
            results += s
            logger.info('{} loaded, {} sequences'.format(proteome, len(s)))
            
        except requests.exceptions.ConnectionError as e:
            logger.error("Error Connecting:", e)
            raise SystemExit(e)
        except requests.exceptions.Timeout as e:
            logger.error("Timeout error:", e)
            raise SystemExit(e)
        except requests.exceptions.TooManyRedirects as e:
            logger.error("Too many redirects:", e)
            raise SystemExit(e)
        except requests.exceptions.RequestException as e:
            logger.error(e)
            raise SystemExit(e)

    with open(seqfile, "w") as output_handle:
        SeqIO.write(results, output_handle, "fasta")
    


# Building database of clusters with cd-hit
@timeit
def build_clusters(paths, logger):
    logger.info('Clustering sequences with cd-hit ...')

    
    import subprocess
    process = subprocess.Popen([cd_hit, '-i', seqfile, '-o', dbfile, '-c', '0.7', '-n', '5'],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        logger.error(stderr)
    logger.debug(stdout)

# Parsing clusters data from .clstr file
def parse_clusters(paths, logger):
    clusters = {}
    with open('{}.clstr'.format(dbfile), 'r') as file:
        for line in file:
            if line[0] == '>':
                name = line[1:-1].replace(' ', '')
                clusters[name] = []
            else:
                seq_name = line.split()[2].split('|')[1]
                clusters[name].append(seq_name)
        logger.info('{} clusters loaded'.format(len(clusters)))
    return clusters

# Building gene name to sequence and gene name to organism dictionaries
def build_dictioaries(results):
    results=SeqIO.parse(results, 'fasta')
    org = {}
    all_taxa = set()
    sequences = {}
    for i in results:
        name = i.name.split('|')[1]
        tax = ''.join(j + '_' for j in i.description.split('OS=')[1].split()[:2])
        org[name] = tax[:-1]
        all_taxa.add(tax[:-1])
        sequences[name] = i
    return sequences, all_taxa, org

# Filter out paralogus
#TODO: more objective method for paralogus filtering, idea - iterative alignment
def filter_paralogus(filtered, org):
    paralogus_free = {} # map clusters to lists of sequences without paralogus
    for key, value in filtered.items():
        used_species= []
        seq_list = []
        for seq in value:
            if org[seq] not in used_species:
                seq_list.append(seq)
                used_species.append(org[seq])
        paralogus_free[key] = seq_list
    return paralogus_free

# Build MSA for paralogus free clusters
@timeit
def msa_one_to_one(paralogus_free, paths, logger):
    logger.info('Building MSA for 1 to 1 clusters')


    for k in paralogus_free:
        cline = ClustalOmegaCommandline(clustal_exec, infile='{}/{}'.format(clusters_dir, k),
                                        outfile='{}/{}'.format(msa_dir, k))
        stdout, stderr = cline()
        if stderr:
            logger.error(stderr)
        logger.debug(stdout)


def write_fatsa(paralogus_free, paths, sequences, org):
    for key, value in paralogus_free.items():
        with open('{}/{}'.format(clusters_dir, key), "w") as file:
            for i in value:
                seq = sequences[i]
                seq.id = org[i]
            SeqIO.write([sequences[x] for x in value], file, "fasta")

def write_fatsa_paralogus(filtered, paths, sequences, org):
    for key, value in filtered.items():
        with open('{}/{}_paralogus'.format(clusters_dir, key), "w") as file:
            for i in value:
                seq = sequences[i]
                seq.id = i + '_' + org[i]
            SeqIO.write([sequences[x] for x in value], file, "fasta")

# Build MSA for clusters with paralogus sequences
@timeit
def msa_paralogus(filtered,paths,logger):
    logger.info('Building MSA for clusters with paralogus')
    for k in filtered:
        cline = ClustalOmegaCommandline(clustal_exec, infile='{}/{}_paralogus'.format(clusters_par_dir, k),
                                        outfile='{}/{}_paralogus'.format(msa_par, k))
        stdout, stderr = cline()
        if stderr:
            logger.error(stderr)
        logger.debug(stdout)

# Build ML tree with IQTREE
@timeit
def ml_trees(paths, paralogus_free, outname, nt, logger):
    for i in paralogus_free:
        logger.info('Building ML tree for {}'.format(i))
        process = subprocess.Popen([iqtree_exec, '-s', '{}/{}'.format(msa_dir, i), '-o', outname,
                                    '-pre', '{}/{}'.format(tree_dir, i), '-nt', nt, '-seed', '12345'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            logger.error(stderr)
        logger.debug(stdout)

# Marge all trees into one file
def marge_trees(paths,paralogus_free):
    with open('{}/MLall.treefile'.format(tree_dir), 'w') as all_tree:
        for i in paralogus_free.keys():
            with open('{}/{}.treefile'.format(tree_dir, i), 'r') as file:
                all_tree.write(file.read())
                all_tree.write('\n')

# Build consensus tree
@timeit
def build_consensus(paths, logger):
    process = subprocess.Popen(
        [iqtree_exec, '-t', '{}/MLall.treefile'.format(tree_dir), '-con', '-pre', '{}/consensus'.format(tree_dir)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        logger.error(stderr)
    logger.debug(stdout)

# Build  MLtrees in IQTREE for clusters with paralogus
@timeit
def ml_trees_para(iqtree_exec, msa_par, filtered, tree_par, nt, logger):
    for i in filtered:
        logger.info('Building ML tree for {} with paralogus'.format(i))
        process = subprocess.Popen([iqtree_exec, '-s', '{}/{}_paralogus'.format(msa_par, i), '-pre',
                                    '{}/{}_paralogus'.format(tree_par, i), '-nt', nt, '-seed', '12345'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            logger.error(stderr)
        logger.debug(stdout)

# Build super tree
@timeit

def marge_paralogus(paths, filtered):
    # Cut genes name from tree with paralogus
    pattern = re.compile('[A-Z0-9]+_')
    with open('{}/Paralogus_all.treefile'.format(tree_par), 'w') as all_tree:
        for k in filtered:
            with open('{}/{}_paralogus.treefile'.format(tree_par, k)) as file:
                all_tree.write(file.read())
                all_tree.write('\n')

    with open('{}/Paralogus_all.treefile'.format(tree_par), 'r') as all_tree:
        with open('{}/all_ng.treefile'.format(tree_par), 'w') as file2:
            trees = all_tree.read()
            m = re.findall(pattern, trees)
            for i in m:
                trees = trees.replace(i, '')
            file2.write(trees)


def build_super_tree(paths, logger):


    # Create .txt file with commands for Clann 
    with open('comands.txt', 'w') as file:
        file.write('hs')
    logger.info('Building super tree tree')
    process = subprocess.Popen([clann_exec, '-n', '-c', 'comands.txt', '{}/all_ng.treefile'.format(tree_par)],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        logger.error(stderr)
    logger.debug(stdout)


# Build ML trees with boostrap in IQTREE
@timeit
def build_boots_tree(paths,outname, filtered, nt, logger):
    for i in filtered:
        logger.info('Building ML tree with bootstrap for {}'.format(i))
        process = subprocess.Popen([iqtree_exec, '-s', '{}/{}'.format(msa_dir, i), '-o', outname, '-bb',
                                    '1000', '-pre', '{}/{}_boots'.format(tree_boots, i), '-nt', nt, '-seed', '12345'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stderr:
            logger.error(stderr)
        logger.debug(stdout)

@timeit
def build_boots_cons(paths,filtered, logger ):
    threshold = 50
    # Extract bootstrap info from trees
    with open('{}/all_boots.treefile'.format(tree_boots), 'w') as all_file:
        with open('{}/filtered_boots.treefile'.format(tree_boots), 'w') as filtered_file:
            for k in filtered:
                with open('{}/{}_boots.contree'.format(tree_boots, k)) as file:
                    f = file.read()
                    all_file.write(f)
                    pattern = re.compile('[0-9]+:')
                    m = re.findall(pattern, f)
                    x = [int(i[:-1]) for i in m]
                    if all([i >= threshold for i in x]):  # Take only trees with support greater than threshold
                        filtered_file.write(f)
                        print('przeszlo')

    # Build consensus for best bootstrap trees
    process = subprocess.Popen(
        [iqtree_exec, '-t', '{}/filtered_boots.treefile'.format(tree_boots), '-con', '-pre',
         '{}/bootscons'.format(tree_boots)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if stderr:
        logger.error(stderr)
    logger.debug(stdout)



def clean_supertree():
    with open('Heuristic_result.txt') as file:
        with open('super_tree_ur', 'w') as file2:
            text = file.read()
            pattern = '\[\d*\.\d*\]'
            text2 = re.sub(pattern, '', text)
            file2.write(text2)


# Write consensus and super trees in Newick format
def write_newick(paths, outname ):
    tree = Phylo.read('{}/consensus.contree'.format(paths.tree_dir), "newick")
    tree.root_with_outgroup({"name": outname})
    Phylo.write(tree, '{}/ML_consensus'.format(paths.base), "newick")

    tree = Phylo.read('{}/bootscons.contree'.format(paths.tree_boots), "newick")
    tree.root_with_outgroup({"name": outname})
    Phylo.write(tree, '{}/bootstrap_consens'.format(paths.base), "newick")

    trees = Phylo.parse('super_tree_ur', "newick")
    tree = next(trees)  # If there is more than one super tree take first one
    tree.root_with_outgroup({"name": outname})
    Phylo.write(tree, '{}/super_tree'.format(paths.base), "newick")


# Clean files
def remove_0_length(tree_file):
    with open(tree_file, 'r') as file:
        text = file.read()
        pattern = ':0.00000'
        text2 = re.sub(pattern, '', text)
    with open(tree_file, 'w') as file:
        file.write(text2)



class PathsHandler():
    
    def __init__(self, args) -> None:
        self.iqtree=args.i
        self.clann=args.s
        self.clustal=args.cl
        self.cd_hit=args.cd
        self.base=args.o
        self.tmp=self.base+'/tmp'
        self.seqfile=self.tmp+"/proteomes.fasta"
        self.dbfile = self.tmp + '/dbout',
        self.tree_dir = self.tmp + '/ml_trees',
        self.msa_dir = self.tmp + '/msa',
        self.clusters_dir = self.tmp + '/clusters',
        self.clusters_par_dir = self.tmp + '/clusters_par',
        self.msa_par = self.tmp + '/msa_par',
        self.tree_par = self.tmp + '/tree_par',
        self.tree_boots = self.tmp + '/tree_boots'




def build_dir_structure(paths, logger):
    try:
        os.mkdir(paths.base)
    except FileExistsError as e:
        logger.error(e)
        exit(1)
    os.mkdir(paths.tmp)
    os.mkdir(paths.clusters_dir)
    os.mkdir(paths.tree_dir)
    os.mkdir(paths.msa_dir)
    os.mkdir(paths.clusters_par_dir)
    os.mkdir(paths.msa_par)
    os.mkdir(paths.tree_par)
    os.mkdir(paths.tree_boots)

def main():
    parser = console_args()
    args = parser.parse_args()
    logger=set_logger()    
    taxa = args.t
    outname = args.og
    nt = str(args.n)
    Paths=PathsHandler(args)
    build_dir_structure(paths, logger)
    download_data(taxa,outname, Paths,logger)
    build_clusters(Paths, logger)
    
    clusters = parse_clusters( Paths,logger)
    sequences, all_taxa, org = build_dictioaries(Paths)
    filtered = {key: value for (key, value) in clusters.items() if set([org[seq] for seq in value]) == all_taxa}
    logger.info(' {} clusters were selected for further analysis'.format(len(filtered)))
    paralogus_free = filter_paralogus(filtered, org)
    write_fatsa(paralogus_free, Paths, sequences, org)
    write_fatsa_paralogus(filtered, Paths, sequences, org)
    msa_one_to_one(paralogus_free, Paths, logger)
    msa_paralogus(filtered, Paths,logger)

    ml_trees(Paths, paralogus_free, outname, nt, logger)
    marge_trees(Paths,paralogus_free)
    build_consensus(Paths, logger)
    ml_trees_para(Paths,filtered,nt, logger)
    marge_paralogus(Paths, filtered)
    build_super_tree(Paths, logger)
    clean_supertree()
    build_boots_tree(Paths, outname, paralogus_free, nt, logger)
    build_boots_cons(Paths,filtered, logger)
    write_newick(Paths, outname )

    # Clear the output files
    remove_0_length('{}/super_tree'.format(Paths))
    remove_0_length('{}/bootstrap_consens'.format(Paths))
    remove_0_length('{}/ML_consensus'.format(Paths))

    # Move trees from tmp file 
    shutil.move('{}/all_boots.treefile'.format(Paths.tree_boots), Paths.base)
    shutil.move('{}/Paralogus_all.treefile'.format(Paths.tree_par), Paths.base)
    shutil.move('{}/MLall.treefile'.format(Paths.tree_dir), Paths.base)
    # Remove temporary files if not in debugg mode
    if not args.D:
        shutil.rmtree(Paths.tmp)    


if __name__== "__main__":
    main()