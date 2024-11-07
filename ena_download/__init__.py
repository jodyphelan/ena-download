"""
This package facilitates the download of data from the ENA in fastq format.
To use it, you need to provide the accession number of the data you want to download.
"""
__version__ = '0.1.3'
import requests
import os
import subprocess as sp
import argparse
from typing import List
import sys
import json
import logging
logging.basicConfig(level=logging.INFO)

def is_valid_accession(accession: str) -> bool:
    """
    Get the URL of the data to download.

    Parameters
    ----------
    accession :
        The run accession number of the data to download.

    Returns
    -------
    str
        The URL of the data to download.

    Examples
    --------
    >>> is_valid_accession("ERR11466368")
    True
    >>> is_valid_accession("ERR0000000")
    Traceback (most recent call last):
    ValueError: Invalid accession number: ERR0000000
    """
    logging.debug(f"Checking if {accession} is a valid accession number")

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    parameters = {
        "result": "read_run",
        "includeAccessions": accession,
        "limit": 10,
        "format": "json"
    }
    response = requests.get(url, params=parameters)
    
    data = json.loads(response.text)

    if len(data) == 0:
        raise ValueError(f"Invalid accession number: {accession}")
    return True

def extract_data_path(accession: str) -> List[str]:
    """
    Get the URL of the data to download.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.

    Returns
    -------
    str
        The URLs of the data to download.

    Examples
    --------
    >>> extract_data_path("ERR11466368")
    ['ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/068/ERR11466368/ERR11466368_1.fastq.gz', 'ftp.sra.ebi.ac.uk/vol1/fastq/ERR114/068/ERR11466368/ERR11466368_2.fastq.gz']
    """
    
    logging.debug(f"Extracting data path for {accession}")

    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    parameters = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
        "format": "json"
    }

    response = requests.get(url, params=parameters)
    if response.status_code != 200:
        raise ValueError(f"Invalid URL: {url}")
    
    logging.debug(f"Response: {response.text}")
    data = json.loads(response.text)
    logging.debug(f"Data found for {accession}: {data}")

    if len(data[0]['fastq_ftp']) == 0:
        raise ValueError(f"No data found for {accession}")

    files = []
    for d in data:
        files += d['fastq_ftp'].split(";")
    
    if len(files) == 0:
        raise ValueError(f"No data found for {accession}")
    
    return files

import signal, os

def handler(signum, frame):
    """Signal handler for the download timeout."""
    raise TimeoutError(f'Download timeout reached, trying again!')

def download_data(accession: str, urls: List[str],timeout: int = 300) -> None:
    """
    Download data from the ENA.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.
    urls : str
        The URLs of the data to download.
    timeout : int
        The timeout in seconds for the download to complete. Default is 300 seconds.

    Returns
    -------
    None
    """

    logging.debug(f"Downloading data for {accession}")
    if not os.path.exists(accession):
        os.mkdir(accession)

    home = os.path.expanduser("~")
    ascp = os.path.join(home, '.aspera/cli/bin/ascp')
    opensshfile = os.path.join(home, '.aspera/cli/etc/asperaweb_id_dsa.openssh')

    for url in urls:
        i=0
        while True:
            sys.stderr.write(f"Attempt {i+1} at downloading {url}...\n")
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            try:
                path = url.replace('ftp.sra.ebi.ac.uk/', 'era-fasp@fasp.sra.ebi.ac.uk:')
                sp.run([
                    ascp, '-T', '-l', '300m', '-P', '33001', '-i', opensshfile, 
                    path, accession + '/'
                ], check=True)
                signal.alarm(0)
                break
            except:
                i+=1
                if i==3:
                    raise TimeoutError(f"Download failed after 3 attempts")
                continue
    
    if accession.startswith('SAM'):
        forward_reads = sorted([f'{accession}/{f}' for f in os.listdir(accession) if f.endswith('_1.fastq.gz')])
        reverse_reads = sorted([f'{accession}/{f}' for f in os.listdir(accession) if f.endswith('_2.fastq.gz')])
        if len(forward_reads) == 0 or len(reverse_reads) == 0:
            raise ValueError(f"Download failed for {accession}")

        sp.run(f"cat {' '.join(forward_reads)} > {os.path.join(accession, accession + '_1.fastq.gz')}", shell=True, check=True)
        sp.run(f"cat {' '.join(reverse_reads)} > {os.path.join(accession, accession + '_2.fastq.gz')}", shell=True, check=True)

    return None

def main(accession: str, timeout: int = 300) -> None:
    """
    Function that calls all the other functions to download data from the ENA.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.
    timeout : int
        The timeout in seconds for the download to complete. Default is 300 seconds.

    Returns
    -------
    None
    """
    
    is_valid_accession(accession)
    
    paths = extract_data_path(accession)

    download_data(accession, paths, timeout)

    
    return None

def cli():
    """
    Entry point for the command line interface. This function is called when the package is called from the command line.
    It uses the argparse package to parse the arguments passed to the command line.

    Returns
    -------
    None
    """
    argparser = argparse.ArgumentParser(description='ENA Download')
    argparser.add_argument('accession', type=str, help='Accession number of the data to download')
    argparser.add_argument('--timeout', default=300, type=int, help='Timeout in seconds for the download to complete. Default is 300 seconds.')
    argparser.add_argument('--debug', action='store_true', help='Print debug information')
    args = argparser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    main(args.accession,args.timeout)
