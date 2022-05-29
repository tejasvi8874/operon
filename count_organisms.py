from helpers import species_list, curl_output, string_id_n_refseq_pairs, stringdb_aliases
from threading import get_ident, Thread
from collections import defaultdict
from queue import SimpleQueue, Empty
import requests
from itertools import chain, tee, groupby
from time import sleep
import sys
from json import loads
import re
from urllib.parse import quote_plus
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

data_path = Path('count_organisms.pkl')
if data_path.exists():
    with open(data_path) as f:
        data = pickle.load(f)
else:
    data = {}
count = 0

def pairwise(iterable):
    # Python 3.10 onwards pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

sessions = defaultdict(lambda: requests.Session())
with open('nofind.txt') as f:
    nofind = {l.split(' ')[-1] for l in f}

def f(genome_organism_id, organism_selection):
    if organism_selection in data:
        return
    try:
        if genome_organism_id in nofind:
            data.setdefault(organism_selection, None)
            return
        try:
            string_refseq_gen = chain(string_id_n_refseq_pairs(genome_organism_id), pairwise(k for k, _ in groupby(m.groups()[0] for m in re.finditer(r"^\d+\.(\S*)\t.*$", stringdb_aliases(genome_organism_id), re.MULTILINE))))
        except:
            print("skip out")
            return
        for _  in range(3):
            resp = sessions[get_ident()].post(
                'https://patricbrc.org/api/genome_feature',
                headers={'Content-Type': 'application/x-www-form-urlencoded'},
                data=f"and(keyword(%22{genome_organism_id}%22),or({','.join(['keyword(%22' + a_string_id + '%22),keyword(%22' + a_refseq + '%22)' for _, (a_string_id, a_refseq) in zip(range(20), string_refseq_gen)])}))&limit(1)"
                )
            resp.raise_for_status()
            features = loads(resp.content)
            if features:
                if type(features) != list:
                    print(features)
                genome_id = features[0]['genome_id']
                data.setdefault(organism_selection, set()).add(genome_id)
                data.setdefault("errors", {}).pop(organism_selection, None)
                global count
                count += 1
                break
        else:
            data.setdefault(organism_selection, None)
            print("Couldn't find", organism_selection, genome_organism_id, file=sys.stderr)
    except Exception as e:
        print("exception", organism_selection, e)
        data.setdefault("errors", {})[organism_selection] = str(e)
        raise

def saver():
    with open(data_path, 'wb') as file:
        pickle.dump(data, file)
def save_trigger():
    while True:
        Thread(target=saver).start()
        sleep(5)
Thread(target=save_trigger, daemon=True).start()

sl = species_list()
sync = 0
if sync:
    for s in sl:
        f(*s)
else:
    with ThreadPoolExecutor() as ex:
        for i, r in enumerate(as_completed([ex.submit(f, *s) for s in sl])):
            r.result()
            print(round((i+1)/len(sl), 1), end='\r')

print("errors", data.get("errors"))
print("Total count", count)
