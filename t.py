from shutil import rmtree
from json import *
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def f(folder):
    try:
        src = folder.joinpath('compare_region')
        if src.exists():
            tar = folder.joinpath('compare_region.json')
            if not tar.exists():
                tar.write_text(dumps({f.stem: loads(f.read_bytes()) for f in src.iterdir()}))
                print("ex")
    except Exception as e:
        print(e)
        raise
#with ProcessPoolExecutor() as exc:
#    exc.map(f, Path('.json_files').iterdir())
for x in Path('.json_files').iterdir():
    p = x.joinpath('compare_region')
    if p.exists():
        rmtree(p)
    #f(x)
