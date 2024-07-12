# maricultureMAGs

__code for calculating mapped-read abundance for MAGs__

-- 1.make_cov_mp.py
```
Usage: generate relative abundance table with multiple cores

Options:
  -h, --help            show this help message and exit
  -i INPUTF, --input=INPUTF
                        '*bam' data path
  -f FASTA, --fasta=FASTA
                        contig fasta file
  -o OUTPUTF, --output=OUTPUTF
                        coverage file name
  -t WORKER, --threads=WORKER
                        working cpu number
  -b BEDT, --bedtools=BEDT
                        bedtools absolute path
```

--2.amr_mapped_range.py
```
Usage: get AMR coverage file from bamfiles

Options:
  -h, --help            show this help message and exit
  -i INPUTF, --input=INPUTF
                        '*bam' data path
  -a AMRFINDER, --amr=AMRFINDER
                        amrfinder result file folder
  -o OUTPUTF, --output=OUTPUTF
                        coverage file name
  -t WORKER, --threads=WORKER
                        working cpu number
  -b BEDT, --bedtools=BEDT
                        bedtools absolute path
```
