#! ${params.python3}

import re
import sys
import requests

url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=${projectID}"
desc = requests.get(url).text

regex_strings = ['WGS|shotgun|RNA-Seq', 'ILLUMINA', 'METAGENOMIC|metagenome|TRANSCRIPTOMIC']
with open('runs.txt', 'w') as out:
    for line in desc.split('\\n'):
        if all(re.search(pat, line) for pat in regex_strings):
            print(line.split(",")[0], file=out)
