import os
import copy
import time
import multiprocessing as mp
import pandas as pd
import optparse
import subprocess
from collections import defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP

def make_cov(path,i,amrD,TotalD,bedtools):
    cmd=[bedtools,"genomecov","-ibam",os.path.join(path,i),"-bg"]
    Ct=defaultdict(int)
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,encoding='utf8')
    while True:
        line = proc.stdout.readline()
        if not line: break
        line=line.rstrip('\n').split()

        # check genomecov file format
        if line[3].isalpha():
            print ('check genomecov results: ',line[3])
            continue
        try:
            dp=int(line[3])
        except:
            dp=float(line[3])

        Ct[line[0]]+=float(int(int(line[2])-int(line[1]))*dp)/TotalL[line[0]]
    TotalD[i.split('.')[0]]=Ct

def main(PATH,Fa,Output,worker,bedtools):
    start_time=time.time()
    TotalD=defaultdict(dict)
    TotalL=defaultdict(int)

    with open(Fa,'r') as F:
        for t,s in SFP(F):
            t=t.split()[0]
            S=len(s)
            TotalL[t]=S
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
            proc=mp.Process(target=make_cov,args=(PATH,i,TotalD,TotalL,bedtools))
            procs.append(proc)
            proc.start()
        for proc in procs:
            proc.join()
    
    D = copy.deepcopy(TotalD)   
    print("--- %s seconds ---" %(time.time() - start_time))
    df=pd.DataFrame.from_dict(D)
    df.fillna(0, inplace=True).round(2)
    df.to_csv('contigAB.'+Output,sep='\t',index_label='contig')

    MAG={}
    MAGl=defaultdict(int)
    for i in TotalL:
        j=i.split('_')
        MAG[i]=j[0]+'_'+j[1]
        MAGl[j[0]+'_'+j[1]]+=TotalL[i]/1000.0

    dfDctg=df.T.mul(TotalL).T
    dfDctg['MAGid']=dfDctg.index.map(MAG)

    #Mapped read by MAG
    dfDmag=dfDctg.groupby('MAGid').sum()
    #RPK by MAG size
    dfDmag=dfDmag.T.div(MAGl).T
    #TSS relative abundance
    dfDmag=100*dfDmag.div(dfDmag.sum())
    dfDmag.round(2).to_csv(Output,sep='\t',index=True)

if __name__=="__main__":
    usage = "generate MAG's relative abundance table with multiple cores"
    parser = optparse.OptionParser(usage)
    parser.add_option("-i","--input",dest="inputf",
            help="\'*bam\' data path",type="string")
    parser.add_option("-f","--fasta",dest="fasta",
            help="contig fasta file",type="string")
    parser.add_option("-o","--output",dest="outputf",
            help="coverage file name",type="string")
    parser.add_option("-t","--threads",dest="worker",
            help="working cpu number",type="int")
    parser.add_option("-b","--bedtools",dest="bedt",
            help="bedtools absolute path",type="int")
    (opts,args)=parser.parse_args()

    if opts.inputf is None or opts.fasta is None or opts.outputf is None:
        parser.print_help()
    else:
        if opts.worker is None:
            opts.worker=1
        if opts.bedt is None:
            opts.bedt='bedtools'
        main(opts.inputf,opts.fasta,opts.outputf,opts.worker.opts.bedt)
