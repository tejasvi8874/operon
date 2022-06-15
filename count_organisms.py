from helpers import species_list, string_id_n_refseq_pairs, get_session, get_genome_id_name
import traceback
from threading import get_ident, Thread
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

for bacteria_not_archaea in (True, False):
    data_path = Path('count_organisms.pkl' if bacteria_not_archaea else 'count_organisms_archaea.pkl')
    if data_path.exists():
        with open(data_path, 'rb') as f:
            full_data = pickle.load(f)
            data = full_data['data']
    else:
        data = {}
        full_data = {'data': data, 'errors': {}}

    def f(genome_organism_id, organism_selection):
        #if data.get(organism_selection):
        if not data.get(organism_selection):
        #if organism_selection in data:
            return
            full_data.setdefault("errors", {}).pop(organism_selection, None)
            return
        try:
            genome_id, patric_genome_name = get_genome_id_name(genome_organism_id)
            if genome_id:
                genome_id_set = data.setdefault(organism_selection, set())
                genome_id_set.add(genome_id)
                data[patric_genome_name] = genome_id_set
            else:
                print(data[organism_selection])
                #data.setdefault(organism_selection, None)
                print("Couldn't find", organism_selection, genome_organism_id, file=sys.stderr)
            full_data.setdefault("errors", {}).pop(organism_selection, None)
        except Exception as e:
            print("exception", organism_selection, genome_organism_id, traceback.format_exc())
            full_data.setdefault("errors", {})[organism_selection] = str(e)

    def saver():
        with open(data_path, 'wb') as file:
            pickle.dump(full_data, file)
    def save_trigger():
        while True:
            Thread(target=saver, daemon=False).start()
            sleep(5)
    Thread(target=save_trigger, daemon=True).start()

    sl = species_list(bacteria_not_archaea)
    print(len(sl), "organisms")
    sync = 0
    if sync:
        for s in sl:
            f(*s)
    else:
        with ThreadPoolExecutor(max_workers=100) as ex:
            for i, r in enumerate(as_completed([ex.submit(f, *s) for s in sl])):
                r.result()
                print(round((i+1)/len(sl), 1), end='\r')
    saver()
    print("errors", full_data["errors"])
    count = len({gid for v in data.values() if isinstance(v, set) for gid in v})
    print("Total count for", "Bacteria" if bacteria_not_archaea else "Archaea", count)
