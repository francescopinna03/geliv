"""Download public sources explicitly, hash bytes, retain source provenance.
No datasets are redistributed with this toolkit. Run from project root.
"""
import argparse, hashlib, json, urllib.request, zipfile
from pathlib import Path
from datetime import datetime, timezone

SOURCES={
 'card': ('https://fmwww.bc.edu/ec-p/data/wooldridge/card.dta','card.dta'),
 'ak-minimal': ('https://economics.mit.edu/sites/default/files/publications/asciiqob.zip','asciiqob.zip'),
 'ak-full': ('https://economics.mit.edu/sites/default/files/inline-files/NEW7080_1.rar','NEW7080_1.rar')
}
def main():
    p=argparse.ArgumentParser(); p.add_argument('dataset',choices=SOURCES)
    args=p.parse_args(); url,name=SOURCES[args.dataset]
    root=Path(__file__).resolve().parents[1]/'data'; root.mkdir(exist_ok=True)
    dst=root/name
    if not dst.exists():
        req=urllib.request.Request(url,headers={'User-Agent':'GELIV-research/0.1'})
        with urllib.request.urlopen(req,timeout=120) as response:
            payload=response.read()
        dst.write_bytes(payload)
    manifest={'file':name,'source':url,'sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),
              'bytes':dst.stat().st_size,'checked_utc':datetime.now(timezone.utc).isoformat()}
    (root/(name+'.source.json')).write_text(json.dumps(manifest,indent=2))
    if dst.suffix=='.zip':
        with zipfile.ZipFile(dst) as z:
            names=[s for s in z.namelist() if Path(s).name.lower()=='asciiqob.txt']
            if len(names)!=1: raise ValueError('Expected one asciiqob.txt; inspect archive')
            (root/'asciiqob.txt').write_bytes(z.read(names[0]))
    if dst.suffix=='.rar':
        print('Extract the original Stata dataset with a RAR utility. See docs/DATA.md.')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__': main()
