#!/usr/bin/env python3
import argparse,csv,gzip,hashlib,re
from pathlib import Path
DNA=set('ACGT'); COMP=dict(zip('ATCG','TAGC'))
def chrom(p): return re.search(r'chr(\d+)',Path(p).name).group(1)
def keyed(xs): return {chrom(x):x for x in xs}
def canon(pos,a,b):
    d=tuple(sorted((a,b))); c=tuple(sorted((COMP[a],COMP[b]))); q=min(d,c); return (int(pos),*q)
def sha(path,decomp=False):
    h=hashlib.sha256(); op=gzip.open if decomp else open
    with op(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def sample_ids(path):
    z=[x.split() for x in Path(path).read_text().splitlines() if x.strip()]
    assert z[0][:2]==['ID_1','ID_2'] and z[1][:2]==['0','0']
    assert all(x[0]==x[1] for x in z[2:]); return [x[0] for x in z[2:]]
def vcf_header(path):
    with gzip.open(path,'rt') as f:
        for x in f:
            if x.startswith('#CHROM'): return x.rstrip().split('\t')[9:]
    raise SystemExit(f'ERROR no #CHROM in {path}')
def vcf_sites(path):
    s=set(); n=0
    with gzip.open(path,'rt') as f:
        for x in f:
            if x.startswith('#'): continue
            q=x.rstrip().split('\t'); a,b=q[3],q[4]
            if len(a)!=1 or len(b)!=1 or a not in DNA or b not in DNA or a==b: raise SystemExit(f'ERROR non-SNP {q[0]}:{q[1]}')
            k=canon(q[1],a,b)
            if k in s: raise SystemExit(f'ERROR duplicate input site {k}')
            s.add(k); n+=1
    return s,n
def scan_haps(path,c,n_samples,sites,expected):
    n=bad=nonbin=dupe=duppos=unsorted=0; seen=set(); prev=None
    with gzip.open(path,'rt') as f:
        for x in f:
            q=x.split()
            if not q: continue
            if len(q)!=5+2*n_samples or q[0]!=c: raise SystemExit(f'ERROR malformed HAPS chr{c}')
            pos=int(q[2]); a,b=q[3],q[4]; k=canon(pos,a,b)
            bad+=k not in sites; dupe+=k in seen; seen.add(k)
            if prev is not None: unsorted+=pos<prev; duppos+=pos==prev
            prev=pos; nonbin+=sum(v not in ('0','1') for v in q[5:]); n+=1
    if n!=expected: raise SystemExit(f'ERROR chr{c}: HAPS {n} != expected {expected}')
    return n,bad,nonbin,dupe,duppos,unsorted,sha(path),sha(path,True)
def nextrec(f):
    for x in f:
        if not x.startswith('#'): return x.rstrip().split('\t')
    return None
def genotype_check(vcf,haps,n_samples):
    checked=missing=mismatch=matched=0
    with gzip.open(vcf,'rt') as vf,gzip.open(haps,'rt') as hf:
        v=nextrec(vf)
        for x in hf:
            h=x.split(); hp=int(h[2]); hk=canon(hp,h[3],h[4])
            while v is not None and canon(v[1],v[3],v[4])!=hk:
                if int(v[1])>hp: raise SystemExit(f'ERROR phased site missing from VCF {h[0]}:{hp}')
                v=nextrec(vf)
            if v is None: raise SystemExit('ERROR VCF ended before HAPS')
            fmt=v[8].split(':'); gi=fmt.index('GT'); sf=v[9:]
            same=(h[3],h[4])==(v[3],v[4]) or (COMP[h[3]],COMP[h[4]])==(v[3],v[4])
            rev=(h[4],h[3])==(v[3],v[4]) or (COMP[h[4]],COMP[h[3]])==(v[3],v[4])
            if not(same or rev) or len(sf)!=n_samples: raise SystemExit(f'ERROR allele/sample mismatch {v[0]}:{v[1]}')
            hb=h[5:]
            for i,s in enumerate(sf):
                p=s.split(':'); gt=p[gi] if gi<len(p) else '.'
                if '.' in gt: missing+=1; continue
                g=re.split(r'[|/]',gt)
                if len(g)!=2 or any(a not in ('0','1') for a in g): missing+=1; continue
                vd=sum(a=='1' for a in g); hd=(hb[2*i]=='1')+(hb[2*i+1]=='1'); hd=hd if same else 2-hd
                checked+=1; mismatch+=vd!=hd
            matched+=1; v=nextrec(vf)
    return checked,missing,mismatch,matched
def logcheck(path,snps,n,ref):
    t=Path(path).read_text(errors='replace')
    core=[re.search(r'Running time:\s*\d+\s*seconds',t),'Main iteration [20/20]' in t,'Normalising graphs' in t,'Solving haplotypes' in t,re.search(rf'\b{snps}\s+SNPs included\b',t),re.search(rf'\b{n}\s+samples\b',t),re.search(rf'\b{ref}\s+reference haplotypes\b',t),('400 H' in t or re.search(r'\b400\s+states per window\b',t))]
    rt=re.search(r'Running time:\s*(\d+)\s*seconds',t)
    extra=[re.search(r'Seed\s*[:=]\s*15052011\b',t,re.I),re.search(r'\b1\s+thread\b',t,re.I),re.search(r'(?:window|windows)[^\n]*\b0\.5\s*Mb\b',t,re.I),re.search(r'(?:effective population|\bNe\b)[^\n]*\b11418\b',t,re.I)]
    return all(core),int(rt.group(1)) if rt else 0,[bool(x) for x in extra]
def main():
    p=argparse.ArgumentParser()
    for a in ['vcfs','haps','samples','logs']: p.add_argument('--'+a,nargs='+',required=True)
    p.add_argument('--check',required=True); p.add_argument('--benchmark-haps',required=True); p.add_argument('--benchmark-sample',required=True)
    p.add_argument('--n',type=int,default=46); p.add_argument('--ref',type=int,default=1006); p.add_argument('--total',type=int,default=4914283)
    p.add_argument('--summary',required=True); p.add_argument('--by-chr',required=True); p.add_argument('--repro',required=True); p.add_argument('--sample-order',required=True); a=p.parse_args()
    V,H,S,L=map(keyed,[a.vcfs,a.haps,a.samples,a.logs]); C={r['chromosome']:r for r in csv.DictReader(open(a.check),delimiter='\t')}; rows=[]; order=None
    totals=dict(rows=0,bad=0,nonbin=0,dupe=0,duppos=0,unsorted=0,checked=0,missing=0,mismatch=0,logfail=0,runtime=0); extras=[]
    for c in map(str,range(1,23)):
        ci=C[c]; ni=int(ci['ct_maf005_snps']); nr=int(ci['retained_for_phasing']); ne=int(ci['shapeit2_actual_excluded']); vo=vcf_header(V[c]); so=sample_ids(S[c])
        if len(so)!=a.n or so!=vo or (order is not None and so!=order): raise SystemExit(f'ERROR sample order/count chr{c}')
        order=so if order is None else order; sites,nv=vcf_sites(V[c])
        if nv!=ni: raise SystemExit(f'ERROR input count chr{c}')
        hs=scan_haps(H[c],c,a.n,sites,nr); gc=genotype_check(V[c],H[c],a.n); lc=logcheck(L[c],nr,a.n,a.ref)
        if nv-hs[0]!=ne or gc[3]!=nr: raise SystemExit(f'ERROR exclusion/site count chr{c}')
        r=[c,nv,ne,nr,hs[0],len(so),2*len(so),*hs[1:6],gc[0],gc[1],gc[2],'PASS' if lc[0] else 'FAIL',lc[1],*['YES' if x else 'NO' for x in lc[2]],hs[6],hs[7],sha(S[c])]; rows.append(r)
        for k,v in zip(['rows','bad','nonbin','dupe','duppos','unsorted'],[hs[0],*hs[1:6]]): totals[k]+=v
        totals['checked']+=gc[0]; totals['missing']+=gc[1]; totals['mismatch']+=gc[2]; totals['logfail']+=not lc[0]; totals['runtime']+=lc[1]; extras.append(lc[2])
    if totals['rows']!=a.total: raise SystemExit('ERROR total phased SNP count')
    prod=next(r for r in rows if r[0]=='20'); bh=sha(a.benchmark_haps,True); bs=sha(a.benchmark_sample); repro=(prod[-2]==bh and prod[-1]==bs)
    fail=totals['bad'] or totals['nonbin'] or totals['dupe'] or totals['duppos'] or totals['unsorted'] or totals['mismatch'] or totals['logfail'] or not repro
    Path(a.summary).parent.mkdir(parents=True,exist_ok=True); hdr=['chromosome','input_snps','expected_excluded','expected_retained','haps_rows','samples','haplotypes','site_identity_mismatches','nonbinary_values','duplicate_exact_sites','duplicate_positions','out_of_order_positions','nonmissing_genotypes_compared','input_missing_skipped','genotype_dosage_mismatches','log_core','runtime_seconds','seed_logged','thread1_logged','window05_logged','Ne11418_logged','haps_sha256','haps_decompressed_sha256','sample_sha256']
    with open(a.by_chr,'w',newline='') as f: w=csv.writer(f,delimiter='\t'); w.writerow(hdr); w.writerows(rows)
    Path(a.sample_order).write_text('\n'.join(order)+'\n')
    with open(a.repro,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerows([['metric','value'],['chromosome',20],['production_haps_decompressed_sha256',prod[-2]],['benchmark_haps_decompressed_sha256',bh],['haps_content_identical','PASS' if prod[-2]==bh else 'FAIL'],['sample_content_identical','PASS' if prod[-1]==bs else 'FAIL'],['reproducibility_status','PASS' if repro else 'FAIL']])
    with open(a.summary,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['metric','value']); vals=[('audit_status','PASS' if not fail else 'FAIL'),('autosomes',22),('samples',a.n),('haplotypes_per_site',2*a.n),('expected_total_snps',a.total),('observed_total_snps',totals['rows']),('site_identity_mismatches',totals['bad']),('nonbinary_values',totals['nonbin']),('duplicate_exact_sites',totals['dupe']),('duplicate_positions',totals['duppos']),('out_of_order_positions',totals['unsorted']),('nonmissing_genotypes_compared_exhaustively',totals['checked']),('input_missing_skipped',totals['missing']),('genotype_dosage_mismatches',totals['mismatch']),('phase_log_core_failures',totals['logfail']),('all_logs_seed_parsed','PASS' if all(x[0] for x in extras) else 'REVIEW'),('all_logs_thread1_parsed','PASS' if all(x[1] for x in extras) else 'REVIEW'),('all_logs_window05_parsed','PASS' if all(x[2] for x in extras) else 'REVIEW'),('all_logs_Ne11418_parsed','PASS' if all(x[3] for x in extras) else 'REVIEW'),('sum_shapeit2_runtime_seconds',totals['runtime']),('chr20_repeatability','PASS' if repro else 'FAIL'),('lassi_structural_readiness','PASS' if not fail else 'FAIL'),('accuracy_limitation','Internal consistency cannot estimate switch error without independent truth or an external-method sensitivity analysis.')]; w.writerows(vals)
    if fail: raise SystemExit('ERROR post-phasing audit failed; inspect outputs')
if __name__=='__main__': main()
