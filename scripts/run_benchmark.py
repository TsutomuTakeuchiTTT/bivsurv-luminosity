#!/usr/bin/env python3
from pathlib import Path
import sys,os,argparse
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from bivsurv.benchmark import run
p=argparse.ArgumentParser();p.add_argument('--config',default='configs/visual_benchmark.json');p.add_argument('--output',default='outputs/visual');p.add_argument('--replicates',type=int);p.add_argument('--sizes',type=int,nargs='+');a=p.parse_args()
print(run(a.config,a.output,a.replicates,a.sizes))
