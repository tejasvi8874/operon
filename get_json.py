import requests
from threading import Thread
from email_validator import validate_email, EmailNotValidError
import sys
import asyncio
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
from helpers import Wait, get_session, get_compare_region_data, get_compare_region_json_path, send_alert

import streamlit as st
from JsonToCoordinates import parse_string_scores


def get_operon_progress_path(genome_id):
    return Path(f'.json_files/{genome_id}/operons_progress')

def get_operon_path(genome_id):
    return Path(f'.json_files/{genome_id}/operons.json')

def operons_in_progress(genome_id):
    try:
        with Wait('.lock_'+genome_id):
            return False
    except PidFileError:
        return True

class LruDict:
    def __init__(self, size):
        self.size = size
        self.d = {}
    def __getitem__(self, key):
        self.d[key] = self.d.pop(key)
        return self.d[key]
    def __setitem__(self, key, val):
        if len(self.d) >= self.size:
            del self.d[next(iter(self.d))]
        self.d[key] = val
    def __contains__(self, key):
        return key in self.d
    def __iter__(self):
        return iter(self.d)


operon_probs_cache = LruDict(128)
def operon_probs(genome_id: str, pegs: frozenset) -> dict[str, float]:
    if genome_id not in operon_probs_cache:
        makedirs(f'.json_files/{genome_id}', exist_ok=True)
        predict_json = get_operon_path(genome_id)
        if not predict_json.exists():
            placeholder = st.empty()
            placeholder.info(f"Please wait while we fetch the data and predict operons. It might take upto {round(len(pegs)/5/60)} minutes.")
            email = placeholder.input("Get email alert on completion", value=st.session.get("email", ""), placeholder='Enter email address', help="Make sure to check spam/junk folder. Email will be recieved from spklab.iitg@gmail.com").strip()
            if email:
                try:
                    email = email_validator(email).email
                except EmailNotValidError as e:
                    placeholder.error(f"Invalid email address\n\n{e}")
            progress_bar = placeholder.progress(0.05)
            progress_file = get_operon_progress_path(genome_id)

            if not operons_in_progress(genome_id):
                get_operons_background(genome_id, pegs)
            while not predict_json.exists() and operons_in_progress(genome_id):
                progress_bar.progress(float(progress_file.read_text()))
                sleep(1)

            if predict_json.exists():
                if email:
                    sent_emails = st.session.setdefault("sent_emails", LruDict(2048))
                    if (email, genome_id) not in sent_emails:
                        sent_emails.add(email, genome_id)
                        Thread(send_alert, email, genome_id).start()
            else:
                st.error(f"Some error occured, please retry and report the genome id to {source_email}")
                raise Exception(f"Error with {genome_id=}")
            placeholder.empty()

        operons = loads(predict_json.read_bytes())
        operon_probs_cache[genome_id] = {int(gene_id): prob for gene_id, prob in operons.items()}
    probs = operon_probs_cache[genome_id]
    return probs

class LoggedThread(Thread):
    def __init__(self, f, *a, **k):
        def wrap_f(*fa, **fk):
            try:
                f(*fa, **fk)
            except:
                err_msg = traceback.format_exc()
                print(err_msg)
                if environ.get('PROD'):
                    send_alert_background(error_email, genome_id, err_msg)
                raise
        super().__init__(wrap_f, *a, **k)

def get_operons_background(genome_id:str, pegs: frozenset) -> dict[str, float]:
    LoggedThread(get_operons, genome_id, pegs).start()

def get_operons(genome_id:str, pegs: frozenset) -> dict[str, float]:
    with Wait('.lock_'+genome_id):
        predict_json = get_operon_path(genome_id)
        if predict_json.exists():
            return loads(predict_json.read_bytes())

        progress_file = get_operon_progress_path(genome_id)
        progress_writer = lambda progress: progress_file.write_text(str(progress))

        progress_writer(0.0)

        compare_region_json_path = get_compare_region_json_path(genome_id)
        if not compare_region_json_path.exists():
            Path(test_operons_path).unlink(missing_ok=True)
        compare_region_data = get_compare_region_data(genome_id, pegs, progress_writer)

        progress_writer(0.50)

        from JsonToCoordinates import to_coordinates

        test_operons_path = f"images_custom/test_operons/{genome_id}"

        if not Path(test_operons_path).exists() or len(list(Path(test_operons_path).glob('*.jpg'))) < len(pegs) - 50:
            coords_filename = to_coordinates(compare_region_data, genome_id)
            print("Coordinates created")

            progress_writer(0.55)
            makedirs(test_operons_path, exist_ok=True)
            run(["java", "CoordsToJpg.java", coords_filename, test_operons_path])
            Path(coords_filename).unlink()
            progress_writer(0.65)


        from test import main
        with Wait('.main_predictor_lock'):
            operons = main(genome_id, progress_writer)

        progress_writer(1.0)

        predict_json.write_bytes(dumps(operons))
        return operons

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
