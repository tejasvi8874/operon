# from pyinstrument import Profiler

# profiler = Profiler()
# profiler.start()

# import cProfile, pstats, io
# from pstats import SortKey
# pr = cProfile.Profile()
# pr.enable()

from gzip import decompress
from io import TextIOWrapper
import numpy as np
from json import dumps, loads
from os import environ
import requests
from urllib.request import urlopen
from pathlib import Path
from contextlib import nullcontext
from subprocess import check_output
import re
from urllib.parse import quote_plus
from typing import Optional
from threading import Thread
from time import time, sleep
from collections import defaultdict

import pandas as pd
from get_json import operon_clusters, operon_probs
import streamlit as st
import sys
import shlex

from helpers import query_keywords, to_pid, curl_output, data, Wait, string_id_n_refseq_pairs, species_list
from pathlib import Path
import shlex
import subprocess
import streamlit.components.v1 as components

if "shell" in st.experimental_get_query_params():

    def run_command(args):
        """Run command, transfer stdout/stderr back into Streamlit and manage error"""
        st.info(f"Running {args}")
        result = subprocess.run(args, capture_output=True, text=True)
        try:
            result.check_returncode()
            st.info(result.stdout)
        except subprocess.CalledProcessError as e:
            st.error(result.stderr)
            raise e
        st.info(f"Finished")
    if st.text_input("Password", type='password') == environ.get("PASSWORD"):
        cmd = st.text_input("Shell", value="ls", placeholder="$ cmd")
        if cmd.startswith("p "):
            st.info(eval(cmd[2:]))
        else:
            run_command(shlex.split(cmd))



tmate_cmd = """bash -ic 'nohup /usr/bin/tmate -S /tmp/tmate.sock new-session -d & disown -a' >/dev/null 2>&1
/usr/bin/tmate -S /tmp/tmate.sock wait tmate-ready
/usr/bin/tmate -S /tmp/tmate.sock display -p "tmate SSH address: #{tmate_ssh}"
/usr/bin/tmate -S /tmp/tmate.sock display -p "tmate web: #{tmate_web}\""""

@st.cache(hash_funcs={TextIOWrapper: lambda _: None})
def setup():
    def data_commit():
        while True:
            time_left = (data.last_update + 60*5) - time()
            if time_left <= 0:
                if data.changed:
                    try:
                        subprocess.check_call(["git", "add", "-A"], cwd=".json_files")
                        try:
                            subprocess.check_call(["git", "commit", "-am", "Update"], cwd=".json_files")
                        except subprocess.CalledProcessError as e:
                            if e.returncode != 1:
                                raise
                        subprocess.check_call(["git", "push"], cwd=".json_files")
                    except subprocess.CalledProcessError as e:
                        print(f"{e.stdout}{e.stderr}", file=sys.stderr)
                        raise
                data.updated(False)
            else:
                sleep(time_left)
    Thread(target=data_commit, name="Git sync").start()

    print("Loading data", file=sys.stderr)
    try:
        for cmd in tmate_cmd.splitlines():
            print(subprocess.check_output(shlex.split(cmd), text=True), file=sys.stderr)
        if not Path('.json_files').exists():
            data_key = "OPERON_DATA_SOURCE"
            if data_key not in environ:
                raise Exception(f"{data_key} environment variable missing. It should contain git repository URL.")
            subprocess.check_call(["git", "clone", "--depth=1", environ["OPERON_DATA_SOURCE"], ".json_files"])
            subprocess.check_call(["git", "config", "--global", "user.email", "operon@git.email"])
            subprocess.check_call(["git", "config", "--global", "user.name", "git.name"])
    except subprocess.CalledProcessError as e:
        print(f"{e.stdout}{e.stderr}", file=sys.stderr)
        raise
streamlit_cloud = environ.get("HOSTNAME", None) == "streamlit"


try:
    st.set_page_config(page_title="Operon Finder", page_icon=":dna:", layout="wide")
except:
    pass
st.write(
    f"<style>{Path('style.css').read_text()}</style>",
    unsafe_allow_html=True,
)

st.title("Operon Finder")

st.markdown("Cluster genes into operons")

st.sidebar.markdown("### Select genome")

manual = "Specify genome ID"
search = "Search genomes"
genome_id_option = st.sidebar.radio("", (search, manual))

if streamlit_cloud:
    with Wait('.setup_lock'):
        setup()

genome_id = None
if genome_id_option == search:
        sample_organisms = defaultdict(lambda: None)
        for p in Path(".json_files").glob("*/genome.json"):
            genome_name_file = p.parent.joinpath('genome_name.txt')
            if not genome_name_file.exists():
                genome_name = loads(p.read_bytes())["docs"][0]["genome_name"]
                genome_name_file.write_text(genome_name)
            else:
                genome_name = genome_name_file.read_text()

            sample_organisms[ genome_name ] = p.parent.name
        for species_name in species_list():
            sample_organisms.setdefault(species_name, None)

        # Prevent model inference on local machine
        # if not streamlit_cloud:
        #     ...
        organism_selection = st.selectbox(
            "Choose organism", sample_organisms, index=0, help="Press Backspace key to change search query"
        )
        genome_id = sample_organisms[organism_selection]

        if not genome_id:
            st.sidebar.error(
                "It may take long to fetch external data for custom organism during first query."
            )

        if not genome_id:
            genome_organism_id = re.search(rb"https://stringdb-static.org/download/protein.links.v11.5/(\d*).protein.links.v11.5.txt.gz", curl_output(f"https://string-db.org/cgi/download?species_text={quote_plus(organism_selection)}")).groups()[0].decode()
            string_refseq_gen = string_id_n_refseq_pairs(genome_organism_id)
            for _, (a_string_id, a_refseq) in zip(range(4), string_refseq_gen):
                features = loads(curl_output('https://patricbrc.org/api/genome_feature' , '--data-raw', f'and(keyword(%22{genome_organism_id}%22),or(keyword(%22{a_string_id}%22),keyword(%22{a_refseq}%22)))&limit(1)'))
                if features:
                    genome_id = features[0]['genome_id']
                    break
            else:
                print(genome_organism_id, a_string_id, a_refseq, file=sys.stderr)
                st.error("No compatible genomes found in PATRIC and STRING database.")
else:
    genome_id = st.sidebar.text_input(
        "Genome ID",
        "262316.17",
        help="Must be available in PATRIC and STRING databases.",
    ).strip()
    if re.match(r"\d+\.\d+", genome_id):
        try:
            for url in (
                f"https://patricbrc.org/api/genome/{genome_id}",
                f"https://stringdb-static.org/download/protein.links.v11.5/{genome_id.split('.')[0]}.protein.links.v11.5.txt.gz",
            ):
                if not requests.head(url).ok:
                    genome_id = None
                    st.sidebar.error(
                        "This genome ID is not supported. Try searching for the organism name instead."
                    )
        except ConnectionError:
            print("Connection error")
    else:
        genome_id = None
        st.sidebar.error("Invalid Genome ID format. E.g. 262316.17")

    
if genome_id:
    # link_button(f"Genome details: {genome_id}", f"https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features")
    st.sidebar.markdown(
        f"**Genome details:** [`{genome_id}`](https://www.patricbrc.org/view/Genome/{genome_id}#view_tab=features)"
    )


st.sidebar.write("---")


s_chk = st.sidebar.checkbox("Run")
submit = s_chk if genome_id else False

def br(times=1):
    for _ in range(times):
        st.write("<br/>", unsafe_allow_html=True)

if genome_id:
    br()

    full_data, sequence_accession_id, gene_locations, approximated_refseqs, refseq_locus_tag_present = to_pid(genome_id)
    if not refseq_locus_tag_present:
        st.warning("The RefSeq approximations used for this genome might reduce the operon prediction accuracy due to absence of the corresponding annotations in the PATRIC database.")
    if not full_data:
        submit = False
        st.sidebar.error(
            "This genome is not supported. Try searching for the organism name instead."
        )
    df = pd.DataFrame.from_dict(
        full_data, orient="index", columns=["RefSeq", "Description", "Protein ID"]
    )
    df.index = df.index.astype(int)
    df = df.sort_index()
    # df.index.rename("PATRIC ID", inplace=True)

    with st.expander("Input genes") if submit else nullcontext():
        if not submit:
            st.markdown("### Input genes")
        st.dataframe(df)

if submit:
    pegs = frozenset(full_data.keys())
    with Wait('.lock_'+genome_id):
        probs = operon_probs(genome_id, pegs)

    operons = []
    with st.expander("Filter operons", True):
        br()
        min_prob = st.slider(
            "Confidence threshold",
            min_value=.0,
            max_value=1.,
            value=0.5,
            step=0.01,
        )

        clusters = operon_clusters(genome_id, pegs, min_prob, probs)
        df["Confidence"] = pd.Series(probs)
        df["Intergenic distance"] = pd.Series([None]*len(df.index))

        # clusters = [{998, 999, 1002}, {1001, 1002, 1003}, {1006, 1007}, {999, 1002, 1010, 1011, 1012}]

        min_len = min((len(c) for c in clusters), default=0)
        max_len = max((len(c) for c in clusters), default=0)

        cluster_size_range = 1, float("inf")
        must_pegs: set[int] = set()
        any_pegs: Optional[set[int]] = None
        keywords: set[str] = set()

        if clusters:
            cluster_size_range = st.slider(
                "Gene count",
                min_value=min_len,
                max_value=max_len,
                value=(min_len, max_len),
                step=1,
            )

            refseq_help = "Comma separated RefSeq IDs"
            refseq_input_label = "Comma separated RefSeq IDs"
            sample_refseqs = df["RefSeq"][:8]
            refseq_prefill = ', '.join(sample_refseqs)
            
            contain_all = st.checkbox("All of the genes",help=refseq_help)
            if contain_all:
                must_pegs_text = st.text_area(
                    refseq_input_label,
                    refseq_prefill,
                    key="all",
                )
                must_pegs = {p.lower() for p in must_pegs_text.split(",")}

            contain_any = st.checkbox("Atleast one of the genes",help=refseq_help)
            if contain_any:
                any_pegs_text = st.text_area(
                    refseq_input_label,
                    refseq_prefill,
                    key="any",
                )
                any_pegs = {p.lower() for p in any_pegs_text.split(",")}

            contain_keyword = st.checkbox("Gene description keywords", value=False, help="Filter operons by contained gene's function descriptions")
            if contain_keyword:
                desc_keyword_txt = st.text_input("", "mce")
                keywords = query_keywords(desc_keyword_txt)

            body: list[str] = []
            for i, cluster in enumerate(clusters):
                if not (
                    cluster_size_range[0] <= len(cluster) <= cluster_size_range[1]
                    and must_pegs.issubset({full_data[j].n_refseq.lower() for j in cluster})
                    and (
                        not any_pegs
                        or any([full_data[j].n_refseq.lower() in any_pegs for j in cluster])
                    )
                    and (
                        not keywords
                        or any(
                            [
                                all(s in full_data[j].desc.lower() for s in keywords)
                                for j in cluster
                            ]
                        )
                    )
                ):
                    continue

                dfx = df.loc[sorted(cluster)]

                # Confidence score is a gene connected to next gene. Last gene will technically have low score so set it to previous score to keep it meaningful
                idx_second_max, idx_max  = np.partition(dfx.index.values, -2)[-2:]
                dfx.loc[idx_max, "Confidence"] = dfx["Confidence"][idx_second_max]

                # TODO: find "next" instead of pid+1 since some pid missing
                for pid in dfx.index:
                    dfx.loc[pid, "Intergenic distance"] = gene_locations[pid+1].start - gene_locations[pid].end if pid+1 in dfx.index else '-'

                operons.append((i, dfx))

    if operons:
        detailed = st.checkbox(f"Detailed view", value=True) 
        st.info(f"{len(operons)} operons found")
        
        save = st.checkbox(f"Save results") 
        if save:
            st.download_button(
                "Download",
                    data="\n".join(
                        ["\t".join(["PATRIC ID", *df.columns.tolist()])]
                        + [
                            f"Operon {i+1}\n" + dfx.to_csv(header=False, sep="\t")
                            for i, dfx in operons
                        ]
                    ),
                file_name=f"{genome_id}-operon.tsv",
            )

        show_all = False
        for i, (operon_num, dfx) in enumerate(operons):
            st.markdown(f"#### Operon {operon_num+1}")
            approximation_note = '"Approximate RefSeq assignment"'
            dfx["RefSeq"] = dfx["RefSeq"].apply(
                    lambda r: f"""<a target="_blank" href="https://www.ncbi.nlm.nih.gov/refseq/?term={r}">{r.upper() + "</a>"
                            + ("<br><a style='text-decoration: none;' target='_self' href='javascript:alert(" + approximation_note + ")'><span title=" + approximation_note + ">⚠️</span></a>" if r in approximated_refseqs else '')
                        }"""
            )
            dfx["Protein ID"] = dfx["Protein ID"].apply(
                lambda r: f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/protein/?term={r}">{r}</a>'
            )
            if not detailed:
                del dfx["Confidence"]
                del dfx["Intergenic distance"]
            st.write(
                dfx.to_html(
                    justify="center",
                    escape=False,
                    classes=["table-borderless"],
                    border=0,
                    formatters={'Confidence': lambda x: f'<b style="background-color: hsl({120*x}, 100%, 75%)">{x:.2f}</b>'} if detailed else None,
                ),
                unsafe_allow_html=True,
            )
            # st.table(dfx)
            if i >= 50 and not show_all:
                show_all = st.checkbox(
                    f"Show remaining {len(operons) - i - 1} operons", value=False
                )
                if not show_all:
                    break

            if detailed:
                c1, _, c2, _ = st.columns([0.05, 0.875, 0.05, 0.025])
                with c1:
                    show_dna = st.button("DNA", key=operon_num)
                if show_dna:
                    start = gene_locations[min(dfx.index)].start-1
                    end = gene_locations[max(dfx.index)].end
                    # https://www.patricbrc.org/view/Genome/511145.12#view_tab=browser&loc=NC_000913%3A63298..63391&tracks=refseqs%2CRefSeqGenes&highlight=
                    fasta = loads(curl_output(
                        f"https://p3.theseed.org/services/data_api/jbrowse/genome/{genome_id}/features/{sequence_accession_id}?reference_sequences_only=false&start={gene_locations[min(dfx.index)].start-1}&end={gene_locations[max(dfx.index)].end}"
                        ))["features"][0]["seq"][:end-start]
                    with c2:
                        st.download_button(label='📥', file_name=f'{genome_id}-operon-{operon_num+1}-dna.txt', key=operon_num, data=fasta)
                    components.html(f"<textarea readonly rows=20 style='width:100%'>{fasta}</textarea>", height=300, scrolling=False)
    else:
        st.error(f"No matching clusters found")

components.html(
    """<script>
/* Components live in their Iframe. Streamlit's context is the first parent.
Find higher parents for each enclosing IFrames.
Streamlit share adds another iframe on top. */
const p = window.parent.parent;
[p, p.parent].forEach(p=>p.postMessage("appLoaded", "*"));</script>
"""
)

# pr.disable()
# s = io.StringIO()
# sortby = SortKey.CUMULATIVE
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print(s.getvalue())

# profiler.stop()
# profiler.print()
