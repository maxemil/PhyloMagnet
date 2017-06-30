#! ${params.python3}

import re
import sys
import requests

url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=${projectID}"
desc = requests.get(url).text

regex_strings = ['WGS|shotgun', 'ILLUMINA', 'METAGENOMIC|metagenome']
with open('runs.txt', 'w') as out:
    for line in desc.split():
        if all(re.search(pat, line) for pat in regex_strings):
            print(line.split(",")[0], file=out)
            #open(line.split(",")[0], 'w').close()
            #print("%s_%s" % ("${projectID}", line.split(",")[0]))
