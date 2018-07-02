#! ${params.python3}
import requests
import shutil
import subprocess


def download_file(url, local_filename):
  r = requests.get(url, stream=True)
  with open(local_filename, 'wb') as f:
      shutil.copyfileobj(r.raw, f)


def get_url_ena(acc):
    url = ""
    if len(acc) == 9:
        url = "http://ftp.sra.ebi.ac.uk/vol1/{}/{}/{}".format(acc[0:3].lower(), acc[0:6], acc)
    elif len(acc) == 10:
        url = "http://ftp.sra.ebi.ac.uk/vol1/{}/{}/{}/{}".format(acc[0:3].lower(), acc[0:6], "00" + acc[-1], acc)
    elif len(acc) == 11:
        url = "http://ftp.sra.ebi.ac.uk/vol1/{}/{}/{}/{}".format(acc[0:3].lower(), acc[0:6], "0" + acc[-2:], acc)
    return url


def get_url_ncbi(acc):
    url = "http://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{}/{}/{}/{}.sra".format(acc[0:3], acc[0:6], acc, acc)
    return url


def main(acc, db):
    url = ""
    if db == "ena":
        url = get_url_ena(acc)
    elif db == "ncbi":
        url = get_url_ncbi(acc)
    local_filename = acc + '.sra'
    download_file(url, local_filename)
    subprocess.call("fastq-dump --gzip --readids --split-spot --skip-technical --clip {}".format(local_filename).split())


if __name__ == '__main__':
    main("$run_id", "${params.database}")
