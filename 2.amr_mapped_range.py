import os
import copy
import time
import multiprocessing as mp
import pandas as pd
import optparse
import subprocess
from glob import glob
from collections import defaultdict

def make_cov(path,i,amrD,TotalD,bedtools):
    cmd=[bedtools,"genomecov","-ibam",os.path.join(path,i),"-d"]
    Cov=defaultdict(int)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf8')
    while True:
        line = proc.stdout.readline()
        if not line: break
        line=line.rstrip('\n').split()
        if line[2].isalpha():
            print ('check genomecov results: ',line[3])
            continue
        dp=int(line[2])
        pos=int(line[1])
        if dp>0:
            if not line[0] in amrD:
                continue
            for j in amrD[line[0]]:
                start=j[0]
                end=j[1]
                gID=j[2]
                gSize=end-start+1
                if start<=pos<=end:
                    Cov[gID]+=100.0*(1.0/gSize)
    TotalD[i.split('.')[0]]=Cov

def get_amr_bed(pathA):
    amrD=defaultdict(list)
    for i in glob(pathA+'/*'):
        ID=i.split('/')[-1]
        with open(i,'r') as A:
            next(A)
            for l in A:
                l=l.rstrip('\n').split('\t')
                amrD[l[1]].append((int(l[2]),int(l[3]),l[0]))
    return amrD

def main(PATH,pathA,Output,worker,bedtools):
    start_time=time.time()
    TotalD=defaultdict(dict)

    amrD=get_amr_bed(pathA)
    dp_F=os.listdir(PATH)
    dp_F=list(filter(lambda x:x.endswith('.bam'),dp_F))
    manager=mp.Manager()
    TotalD=manager.dict()
    chunk=len(dp_F)//int(worker)
    if worker > mp.cpu_count():
        print ('The number of workers has exceeded the number of maximum cores...\n so #workers is changed to #maximum_cores\n')
        worker=mp.cpu_count()
    for n,i in enumerate(range(0,len(dp_F),int(worker))):
        if n==len(dp_F)//int(worker)+1:
            works=dp_F[i:]
        else:
            works=dp_F[i:i+int(worker)]
        procs=[]
        print (works)
        for i in works:
            proc=mp.Process(target=make_cov,args=(PATH,i,amrD,TotalD,bedtools))
            procs.append(proc)
            proc.start()
    
        for proc in procs:
            proc.join()
    
    D = copy.deepcopy(TotalD)   
    print("--- %s seconds ---" %(time.time() - start_time))
    df=pd.DataFrame.from_dict(D)
    df.fillna(0, inplace=True).round(2)
    df.to_csv(Output,sep='\t',index_label='contig')

if __name__=="__main__":
    usage = "get AMR coverage file from bamfiles"
    parser = optparse.OptionParser(usage)
    parser.add_option("-i","--input",dest="inputf",
            help="\'*bam\' data path",type="string")
    parser.add_option("-a","--amr",dest="amrfinder",
            help="amrfinder result file folder",type="string")
    parser.add_option("-o","--output",dest="outputf",
            help="coverage file name",type="string")
    parser.add_option("-t","--threads",dest="worker",
            help="working cpu number",type="int")
    parser.add_option("-b","--bedtools",dest="bedt",
            help="bedtools absolute path",type="int")
    (opts,args)=parser.parse_args()

    if opts.inputf is None or opts.amrfinder is None or opts.outputf is None:
        parser.print_help()
    else:
        if opts.worker is None:
            opts.worker=1
        if opts.bedt is None:
            opts.bedt='bedtools'
        main(opts.inputf,opts.amrfinder,opts.outputf,opts.worker,opts.bedt)
