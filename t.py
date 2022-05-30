from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from json import loads, dumps

root = Path('.json_files')
def f(fol):
    print(fol)
    j = fol.joinpath('genome.json')
    if j.exists():
        print("in")
        d = loads(j.read_bytes())
        j.write_text(dumps(d["docs"]))
        print("out")
with ProcessPoolExecutor() as e:
    e.map(f, root.iterdir())
#for fol in root.iterdir():
#    f(fol)
