from pathlib import Path
from gzip import compress, decompress
from concurrent.futures import ProcessPoolExecutor
from json import loads, dumps

root = Path('.json_files')
def f(fol):
    print(fol)
    j = fol.joinpath('compare_region.json.gz')
    if j.exists():
        print("in")
        d = loads(decompress(j.read_bytes()))
        print("loaded")
        if isinstance(d, dict):
            j.write_bytes(compress(dumps(list(d.values())).encode()))
        print("out")
if 1:
    for fol in root.iterdir():
        f(fol)
else:
    with ProcessPoolExecutor() as e:
        e.map(f, root.iterdir())
