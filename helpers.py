from heapq import heappush
from string import ascii_letters
import traceback
from shutil import rmtree
from gzip import compress
from concurrent.futures import ThreadPoolExecutor, as_completed
import pickle
from itertools import chain, groupby, tee
from threading import get_ident, Thread
from multiprocessing import Process
import threading
from typing import Optional, NamedTuple, Iterator
from string import ascii_letters
import sys
from collections import namedtuple, defaultdict
import requests
import re
from time import sleep
from os import makedirs, mkdir
from pathlib import Path
from gzip import decompress
from urllib.request import urlretrieve
from glob import glob
from functools import lru_cache
from subprocess import run, check_output, DEVNULL
from json import dump, dumps, load, loads, JSONDecodeError
from time import time
import smtplib
import codecs
from os import environ
from email.message import EmailMessage


import logging
from logging.handlers import RotatingFileHandler
import time

def get_logger():
    #Setup logger
    logger = logging.getLogger(__name__)
    FORMAT = "[%(asctime)s %(threadName)s %(filename)s->%(funcName)s():%(lineno)s]%(levelname)s: %(message)s"
    logging.basicConfig(format=FORMAT, level=logging.INFO)
    #Log to file
    logging_filename = Path('logs/main.log')
    logging_filename.parent.mkdir(parents=True, exist_ok=True)
    handler = RotatingFileHandler(logging_filename, maxBytes=1000000, backupCount=10) #10 files of 1MB each
    handler.setFormatter(logging.Formatter(FORMAT))
    logger.addHandler(handler)
    return logger
logger = get_logger()


class ServerBusy(Exception):
    pass


@lru_cache(128)
def get_output(url)->bytes:
    resp = get_session().get(url)
    resp.raise_for_status()
    return resp.content


PatricMeta = namedtuple('PatricMeta', ['n_refseq', 'desc', 'protein_id'])
LocInfo = namedtuple('LocInfo', ['start', 'end'])

class PidData(NamedTuple):
    full_data: tuple[dict[int]]
    sequence_accession_id: str
    gene_locations: dict[str, LocInfo]
    approximated_refseqs: Optional[list[str]]
    refseq_locus_tag_present: bool

@lru_cache(128)
def to_pid( genome_id: str) -> PidData:
    approximated_refseqs = []

    genome_organism_id = genome_id.split('.')[0]
    feature_data = get_genome_data(genome_id)

    gene_locations = {}
    full_data = {}
    assert feature_data
    refseq_locus_tag_present = False
    used_stripped_numeric_refseqs = {normalize_refseq(feature.get("refseq_locus_tag") or feature.get("gene", "None")).rstrip(ascii_letters) for feature in feature_data}
    used_stripped_numeric_refseqs.discard("") # If no numbers in there. E.g. https://www.patricbrc.org/view/Genome/214092.191#view_tab=features
    i = 0
    max_i = 2*len(feature_data)
    while i < len(feature_data) and i < max_i:
        feature = feature_data[i]
        i += 1

        patric_id = int(feature["patric_id"].split(".")[-1])
        refseq = feature.get("refseq_locus_tag") or feature.get("gene")
        refseq_locus_tag_present = refseq_locus_tag_present or "refseq_locus_tag" in feature
        protein_id = feature.get("protein_id", "None")

        if refseq:
            n_refseq = normalize_refseq(refseq)
        else: # Some genes have neither refseq locus id nor gene symbol assigned https://patricbrc.org/view/Feature/PATRIC.83332.12.NC_000962.CDS.9963.10160.fwd#view_tab=overview
            for delta in (1, -1):
                if patric_id+delta in full_data:
                    adjacent_full_data = full_data[patric_id+delta]
                    adjacent_num_suffixed_refseq = adjacent_full_data.n_refseq.rstrip(ascii_letters)
                    if not adjacent_num_suffixed_refseq:
                        continue # If no numbers in there. E.g. https://www.patricbrc.org/view/Genome/214092.191#view_tab=features
                    refseq_prefix, adjacent_refseq_counter = get_prefix_counter(adjacent_num_suffixed_refseq)
                    n_refseq = refseq_prefix + str(adjacent_refseq_counter - delta)
                    if n_refseq in used_stripped_numeric_refseqs: # Adjacent refseqs already assigned
                        continue

                    if protein_id == "None":
                        if adjacent_full_data.protein_id[-1].isdigit():
                            protein_prefix, adjacent_protein_counter = get_prefix_counter(adjacent_full_data.protein_id)
                            protein_id = protein_prefix + str(adjacent_protein_counter - delta)
                    approximated_refseqs.append(n_refseq)
                    break
            else:
                feature_data.append(feature)
                continue

        desc: str = feature["product"]
        desc_col_loc = desc.find(': ')
        if desc_col_loc != -1:
            desc = desc[desc_col_loc + 2:]

        full_data[patric_id] = PatricMeta(
            n_refseq=n_refseq, desc=desc, protein_id=protein_id
        )

        gene_locations[patric_id] = LocInfo(start=feature['start'], end=feature['end'])

    # Different genes can have different accession ID (though rare) E.g. first and last gene of 798128.4
    # Prioritize ID present in first gene.
    sequence_accession_id = feature_data[0]["sequence_id"]
    return PidData(full_data, sequence_accession_id, gene_locations, approximated_refseqs, refseq_locus_tag_present)

def query_keywords(query: str) -> set[str]:
    return {qs.lower() for qs in query.split(' ') if qs}


@lru_cache(100)
def get_genome_data(genome_id: str) -> list[dict]:
    genome_data_dir = f'.json_files/{genome_id}'
    genome_data_path = Path(f'{genome_data_dir}/genome.json')
    if genome_data_path.exists():
        genome_data = loads(genome_data_path.read_bytes())
    else:
        #E.g. curl 'https://patricbrc.org/api/genome_feature/?&http_accept=application/solr+json' -H 'Content-Type: application/x-www-form-urlencoded' --data-raw 'rql=eq%28genome_id%252C83332.12%29%2526and%28eq%28feature_type%252C%252522CDS%252522%29%252Ceq%28annotation%252C%252522PATRIC%252522%29%29%2526sort%28%252Bfeature_id%29%2526limit%2825000%29' --compressed
        genome_data = get_session().post('https://patricbrc.org/api/genome_feature/',
            data={ 'rql': f'eq(genome_id%2C{genome_id})%26and(eq(feature_type%2C%2522CDS%2522)%2Ceq(annotation%2C%2522PATRIC%2522))%26sort(%2Bfeature_id)%26limit(25000)' }
        ).json()
        if genome_data:
            makedirs(genome_data_dir, exist_ok=True)
            with open(genome_data_path, 'w') as f:
                dump(genome_data, f)
    return genome_data

sessions = defaultdict(lambda: requests.Session())
session_lock = threading.Lock()
def get_session():
    with session_lock:
        return sessions[get_ident()]

def stringdb_aliases(genome_organism_id) -> str:
    path = Path(f'.json_files/alias/{genome_organism_id}.txt')
    if path.exists():
        return path.read_text()
    aliases = decompress(get_output(f"https://stringdb-static.org/download/protein.aliases.v11.5/{genome_organism_id}.protein.aliases.v11.5.txt.gz")).decode()
    #TODO: Handle saving?
    #path.write_text(aliases)
    return aliases

def string_id_n_refseq_pairs(genome_organism_id: str) -> tuple[str,str]:
    for match in re.finditer(r"^\d+\.(\S*)\t(?:\S*:(\S*)\tBLAST_KEGG_KEGGID|(\S*)\tBLAST_UniProt_GN_(?:OrderedLocusNames|ORFNames))$", stringdb_aliases(genome_organism_id), re.MULTILINE):
        string_id, refseq1, refseq2 = match.groups()
        yield string_id.removeprefix('gene:'), normalize_refseq(refseq1 or refseq2)

normalize_refseq = str.lower

def get_prefix_counter(string):
        assert string and string[-1].isdigit(), f"num suffix required. Got {repr(string)}"
        digits = []
        dot_seen = False
        for c in reversed(string):
                if c.isdigit():
                    digits.append(c)
                elif c == '.' and not dot_seen:
                    digits.append(c)
                    dot_seen = True
                else:
                    break
        while digits and digits[-1] == '0':
            digits.pop()
        if not digits:
            raise Exception(f"Can't get prefix of {string}")
        return string[:-len(digits)], (float if dot_seen else int)(''.join(reversed(digits)))

def species_list() -> list[tuple[str, str]]:
    species_list_path = Path(".json_files/species.json")
    if species_list_path.is_file():
        return loads(species_list_path.read_bytes())
    # Considering only Bacteria for now. Archaea might work too.
    species = sorted(re.findall(r"^(\d+)\t\S+\t[^\t]+\t([^\t]+)\tBacteria$", get_output("https://stringdb-static.org/download/species.v11.5.txt").decode(), re.MULTILINE))
    species_list_path.write_text(dumps(species))
    return species

def pairwise(iterable):
    # Python 3.10 onwards pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def get_genome_id(genome_organism_id) -> Optional[str]:
    string_refseq_gen = chain(string_id_n_refseq_pairs(genome_organism_id), pairwise(k.removeprefix('gene:') for k, _ in groupby(m.groups()[0] for m in re.finditer(r"^\d+\.(\S*)\t.*$", stringdb_aliases(genome_organism_id), re.MULTILINE))))
    for _  in range(3):
        # curl https://patricbrc.org/api/genome_feature --data-raw 'and(keyword(%2283332%22),keyword(%22gene%22))'
        features = get_session().post(
            'https://patricbrc.org/api/genome_feature',
            headers={'Content-Type': 'application/x-www-form-urlencoded'},
            data=f"and(keyword(%22{genome_organism_id}%22),or({','.join(['keyword(%22' + a_string_id + '%22),keyword(%22' + a_refseq + '%22)' for _, (a_string_id, a_refseq) in zip(range(20), string_refseq_gen)])}))&limit(1)"
            ).json()
        if features:
            return features[0]["genome_id"]

def valid_organisms() -> Iterator[tuple[str, Optional[set[str]]]]:
    return [(name, genome_ids) for name, genome_ids in pickle.loads(Path('count_organisms.pkl').read_bytes())['data'].items() if genome_ids]


def get_compare_region_json_path(genome_id):
    return Path(f".json_files/{genome_id}/compare_region.json.gz")

def track_call(orig_func):
    def track_call_wrapper(*a, **k):
        try:
            return orig_func(*a, **k)
        except Exception as e:
            logger.error(e, f"{orig_func.__name__}({a}; {k})")
            raise
    return track_call_wrapper


def get_compare_region_data(genome_id, pegs, progress_clb=None):
    compare_region_json_path = get_compare_region_json_path(genome_id)
    if compare_region_json_path.exists():
        return loads(decompress(compare_region_json_path.read_bytes()))

    gene_figure_name = {f"fig|{genome_id}.peg.{i}" for i in pegs}
    compare_region_temp = compare_region_json_path.parent.joinpath('compare_region')
    compare_region_temp.mkdir(parents=True, exist_ok=True)


    @track_call
    def get_compare_region(fig_gene):
        temp_json_path = compare_region_temp.joinpath(f'{fig_gene}.json')
        if temp_json_path.exists():
            try:
                return loads(temp_json_path.read_bytes())
            except JSONDecodeError as e:
                logger.error('JSONDecodeError', fig_gene, e, temp_json_path.read_bytes(), 'Retrying', fig_gene)
                temp_json_path.unlink()
        data = '{"method": "SEED.compare_regions_for_peg", "params": ["' + fig_gene + '", 5000, 20, "pgfam", "representative+reference"], "id": 1}'
        for _ in range(3):
            resp = get_session().post('https://p3.theseed.org/services/compare_region', data=data)
            if resp.ok:
                break
            msg = resp.content.decode()
            err_msg = f"Error with {fig_gene}. " + """Doing with requests though curl command should be: curl --fail --max-time 300 --data-binary '{"method": "SEED.compare_regions_for_peg", "params": ["{fig_gene}", 5000, 20, "pgfam", "representative+reference"], "id": 1}' https://p3.theseed.org/services/compare_region --compressed\n""" + msg
            logger.error(err_msg)
            if '502 Bad Gateway' in msg:
                sleep(5)
                continue
            resp.raise_for_status()
        else:
            resp.raise_for_status()

        temp_json_path.write_bytes(resp.content)
        return resp.json()

    compare_region_data = []
    with ThreadPoolExecutor(max_workers=25) as executor:
        for i, r in enumerate(as_completed([executor.submit(get_compare_region, g) for g in gene_figure_name])):
            if r.exception():
                raise r.exception()
            compare_region_data.append(r.result())
            progress_clb and progress_clb((i+1)/len(gene_figure_name)*0.50)

    assert len(compare_region_data) == len(gene_figure_name), "Error in compare region data fetch"
    compare_region_json_path.write_bytes(compress(dumps(compare_region_data).encode()))
    rmtree(compare_region_temp)

    return compare_region_data

def logged_background(*, is_process, target, args):
    def wrap_target(*a):
        try:
            target(*a)
        except:
            err_msg = traceback.format_exc()
            logger.critical(err_msg)
            if environ.get('PROD'):
                send_alert_background(error_email, genome_id, err_msg)
            raise
    (Process if is_process else Thread)(target=wrap_target, args=args).start()

source_email = codecs.encode('fcxyno.vvgt@tznvy.pbz', 'rot_13')
error_email = codecs.encode('gfgbzne@bhgybbx.pbz', 'rot_13')

def send_alert_background(dest_email, genome_id, err_msg):
    logged_background(is_process=False, target=send_alert, args=(dest_email, genome_id, err_msg))

def send_alert(dest_email, genome_id, err_msg):
    logger.info("Attempt email Send")

    msg = EmailMessage()
    msg['From'] = source_email
    msg['Subject'] = 'Operon Finder task completed'
    msg['To'] = ', '.join([dest_email,])

    msg.set_content(f'Error with {genome_id}' if err_msg else f'The operon predictions for the requested genome id, {genome_id} are now available <a src="https://apps.streamlitusercontent.com/tejasvi8874/operon/main/web.py/+/?genome_id={83332.12}">here</a>.\n\nRegards,\nSCBL - IIT Guwahati', subtype='html')

    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as smtp:
        smtp.login(source_email, environ.get('EMAIL_PASSWORD'))
        smtp.send_message(msg)

    logger.info("Sent!")

