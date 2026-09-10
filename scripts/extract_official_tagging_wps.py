#!/usr/bin/env python3
"""Extract WP values only; keep the generated configuration out of public Git."""
import argparse, gzip, hashlib, json
from pathlib import Path

def load(path):
    data=path.read_bytes()
    return json.loads(gzip.decompress(data) if path.suffix=='.gz' else data)

def working_points(path):
    return next(c['data'] for c in load(path)['corrections'] if c['name']=='UParTAK4_wp_values')

def value(node, **choices):
    while isinstance(node,dict):
        if node['nodetype']!='category': raise ValueError('Unexpected WP schema')
        node=next(v['value'] for v in node['content'] if v['key']==choices[node['input']])
    return float(node)

def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--btag',type=Path,default=Path('private/btagging.json.gz'))
    p.add_argument('--ctag',type=Path,default=Path('private/ctagging.json.gz'))
    p.add_argument('--output',type=Path,default=Path('data/Tagging/official_2024.txt'))
    a=p.parse_args();b,c=working_points(a.btag),working_points(a.ctag)
    rows=['# UParTAK4 Summer24 NanoAODv15; WP values ONLY, no SF applied.']
    for wp in ['M','T']:
        rows.append(f'b{wp} {value(b,working_point=wp):.17g}')
        for axis in ['CvB','CvL']:
            rows.append(f'c{wp}{axis} {value(c,working_point=wp,axis=axis):.17g}')
    rows += ['# btag_sha256 '+hashlib.sha256(a.btag.read_bytes()).hexdigest(),
             '# ctag_sha256 '+hashlib.sha256(a.ctag.read_bytes()).hexdigest()]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text('\n'.join(rows)+'\n')
    print('Wrote private-derived WP configuration:',a.output)
if __name__=='__main__':main()
