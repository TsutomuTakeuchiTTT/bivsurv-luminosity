#!/usr/bin/env python3
from pathlib import Path
import sys,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bivsurv.plotting import make_all
p=argparse.ArgumentParser();p.add_argument('--output',default='outputs/visual');a=p.parse_args()
print('Generated',len(make_all(a.output)),'independent charts (PDF and PNG).')
