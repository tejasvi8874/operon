import requests
from shutil import rmtree
import sys
import asyncio
import aiohttp
from functools import lru_cache
from gzip import decompress, compress
from shutil import move, rmtree
from os import makedirs
from pathlib import Path
from typing import Optional
from urllib.request import urlretrieve
from glob import glob
from subprocess import run, check_output
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from json import loads, dump, dumps
from helpers import data, Wait, get_session

import streamlit as st
from JsonToCoordinates import parse_string_scores


def get_operons(genome_id:str, pegs: frozenset) -> dict[str, float]:
    placeholder = st.empty()
    placeholder.info("Please wait while we fetch the data and predict operons. It might take upto 15 minutes.")

    progress_bar = st.progress(0.05)
    genome_data_changed = False
    gene_figure_name = {f"fig|{genome_id}.peg.{i}" for i in pegs}

    compare_region_json_path = Path(f".json_files/{genome_id}/compare_region.json.gz")

    if compare_region_json_path.exists():
        compare_region_data = loads(decompress(compare_region_json_path.read_bytes()))
    else:
        genome_data_changed = True
        compare_region_temp = compare_region_json_path.parent.joinpath('compare_region')
        compare_region_temp.mkdir(parents=True, exist_ok=True)

        compare_region_data = []

        compare_region_lock = threading.Lock()
        def get_compare_region(fig_gene):
            temp_json_path = compare_region_temp.joinpath(f'{fig_gene}.json')
            if temp_json_path.exists():
                compare_region_data.append(loads(temp_json_path.read_bytes()))
                return
            data = '{"method": "SEED.compare_regions_for_peg", "params": ["' + fig_gene + '", 5000, 20, "pgfam", "representative+reference"], "id": 1}'
            resp = get_session().post('https://p3.theseed.org/services/compare_region', data=data)
            if not resp.ok:
                err_msg = f"Error with {fig_gene}. " + """Doing with requests though curl command should be: curl --fail --max-time 300 --data-binary '{"method": "SEED.compare_regions_for_peg", "params": ["{fig_gene}", 5000, 20, "pgfam", "representative+reference"], "id": 1}' https://p3.theseed.org/services/compare_region --compressed\n""" + resp.content.decode()
                print(err_msg, file=sys.stderr)
                resp.raise_for_status()
            temp_json_path.write_bytes(resp.content)
            with compare_region_lock:
                compare_region_data.append(resp.json())

        with ThreadPoolExecutor(max_workers=100) as executor:
            for i, r in enumerate(as_completed([executor.submit(get_compare_region, g) for g in gene_figure_name])):
                r.result()
                progress_bar.progress((i+1)/len(gene_figure_name)*0.50)

        assert len(compare_region_data) == len(gene_figure_name), "Error in compare region data fetch"
        compare_region_json_path.write_bytes(compress(dumps(compare_region_data).encode()))
        rmtree(compare_region_temp)

    progress_bar.progress(0.50)

    from JsonToCoordinates import to_coordinates

    test_operons_path = f"images_custom/test_operons/{genome_id}"
    if genome_data_changed:
        data.updated()
        Path(test_operons_path).unlink(missing_ok=True)

    if not Path(test_operons_path).exists() or len(list(Path(test_operons_path).glob('*.jpg'))) < len(pegs) - 50:
        coords_filename = to_coordinates(compare_region_data, genome_id)
        print("Coordinates created")

        progress_bar.progress(0.55)
        makedirs(test_operons_path, exist_ok=True)
        run(["java", "CoordsToJpg.java", coords_filename, test_operons_path])
        Path(coords_filename).unlink()
        progress_bar.progress(0.65)

    placeholder.empty()

    from test import main
    with Wait('.main_predictor_lock'):
        return main(genome_id, progress_bar)

@lru_cache(128)
def operon_probs(genome_id: str, pegs: frozenset) -> dict[str, float]:
    makedirs(f'.json_files/{genome_id}', exist_ok=True)
    predict_json = Path(f'.json_files/{genome_id}/operons.json')
    if predict_json.exists():
        operons = loads(predict_json.read_bytes())
        1+1
    else:
        operons = get_operons(genome_id, pegs)
        with open(predict_json, 'w') as f:
            dump(operons, f)
        data.updated()
    # JSON keys can only be strings
    return {int(gene_id): prob for gene_id, prob in operons.items()}

def operon_clusters(genome_id: str, pegs: frozenset[int], min_prob: float, probs: dict[int, float]) -> list[set[int]]:
    peg_next  = {}
    prev = -1
    for peg in sorted(pegs):
        peg_next[prev] = peg
        prev = peg
    peg_nums = sorted(
        [[peg_num, peg_next[peg_num]] for peg_num, prob in probs.items() if peg_num in peg_next and prob >= min_prob]
    )
    # pair = [op, op[:op.rfind('.')+1] + str(peg_num + 1)]

    clusters: list[set[int]] = []
    for peg_num, next_peg_num in peg_nums:
        if clusters and peg_num in clusters[-1]:
            clusters[-1].add(next_peg_num)
        else:
            clusters.append({peg_num, next_peg_num})

    return clusters


