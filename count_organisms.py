from helpers import species_list, curl_output, string_id_n_refseq_pairs, stringdb_aliases
from itertools import chain, tee
from time import sleep
import sys
from json import loads
import re
from urllib.parse import quote_plus
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

data = {}
count = 0

def pairwise(iterable):
    # pairwise('ABCDEFG') --> AB BC CD DE EF FG
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def f(species_name, genome_organism_id):
    if species_name in data:
        return
    try:
        organism_selection = species_name
        print(f"[{genome_organism_id}]")
        string_refseq_gen = chain(string_id_n_refseq_pairs(genome_organism_id), pairwise(m.groups()[0] for m in re.finditer(r"^\d+\.(\S*)\t.*$", stringdb_aliases(genome_organism_id), re.MULTILINE)))
        for _  in range(3):
            features = loads(curl_output(
                'https://patricbrc.org/api/genome_feature',
                '--data-raw',
                f"and(keyword(%22{genome_organism_id}%22),or({','.join(['keyword(%22' + a_string_id + '%22),keyword(%22' + a_refseq + '%22)' for _, (a_string_id, a_refseq) in zip(range(20), string_refseq_gen)])}))&limit(1)"
                ))
            if features:
                if type(features) != list:
                    print(features)
                genome_id = features[0]['genome_id']
                data.get(species_name, set()).add(genome_id)
                global count
                count += 1
                print(count)
                break
        else:
            print("Couldn't find", organism_selection, genome_organism_id, file=sys.stderr)
    except Exception as e:
        print("exception", species_name, e)
        data.get("errors", {})[species_name] = str(e)
        raise

for s in species_list():
    f(s)
#with ThreadPoolExecutor() as ex:
#    sl = species_list()
#    for i, r in enumerate(as_completed([ex.submit(s) for s in sl])):
#        r.result()
#        print(round((i+1)/len(sl), 1), end='\r')
#
#with open('count_organisms.pkl', 'wb') as file:
#    pickle.dump(data, file)
#print("errors", data.get("errors"))
#print("Total count", count)
