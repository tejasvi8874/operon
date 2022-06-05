from get_json import get_compare_region_data
from helpers import valid_organisms
from concurrent.futures import ThreadPoolExecutor, as_completed

def log(*a):
    with open('compare_region.log', 'a') as f:
        print(*a, file=f)

def g(genome_id):
    try:
        log('Doing', genome_id)
        pegs = [int(feature["patric_id"].split(".")[-1]) for feature in get_genome_data(genome_id)]
        get_compare_region_data(pegs)
    except Exception as e:
        log(genome_id, e)


genome_ids = [genome_id for _, cur_genome_ids in reversed(valid_organisms()) for genome_id in cur_genome_ids]

with ThreadPoolExecutor(max_workers=5) as ex:
    for i, r in enumerate(as_completed([ex.submit(g, genome_id) for genome_id in genome_ids]))
        r.result()
        print(round((i+1)/len(genome_ids), 1), end='\r')
