from helpers import valid_organisms, get_genome_data, get_compare_region_data
from concurrent.futures import ThreadPoolExecutor, as_completed

def g(genome_id):
    try:
        print('Doing', genome_id)
        pegs = [int(feature["patric_id"].split(".")[-1]) for feature in get_genome_data(genome_id)]
        get_compare_region_data(genome_id, pegs)
    except Exception as e:
        print(genome_id, e)
        raise


genome_ids = [genome_id for _, cur_genome_ids in reversed(valid_organisms()) for genome_id in cur_genome_ids]

for genome_id in genome_ids:
    g(genome_id)
