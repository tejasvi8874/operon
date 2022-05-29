"""
from pathlib import Path
root = Path('.json_files/alias')
rem = 0
with open('nofind.txt') as f:
    for oid in [l.split(' ')[-1].strip() for l in f]:
        f = root.joinpath(oid+'.txt')
        if f.exists():
            f.unlink()
            rem += 1
            print(rem)
"""

from pathlib import Path
from gzip import compress
from concurrent.futures import ProcessPoolExecutor

root = Path('.json_files')

def f(fol):
    if '.' in fol.name and fol.is_dir():
        cr = fol.joinpath('compare_region.json')
        tar = fol.joinpath(cr.name + '.gz')
        if not tar.exists():
            tar.write_bytes(compress(cr.read_bytes()))
            cr.unlink()
for fol in root.iterdir():
    f(fol)
#with ProcessPoolExecutor() as ex:
#    for fol in root.iterdir():
#        ex.submit(f, fol)
