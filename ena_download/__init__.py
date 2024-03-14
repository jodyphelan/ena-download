"""
This package facilitates the download of data from the ENA in fastq format.
To use it, you need to provide the accession number of the data you want to download.
"""
__version__ = '0.1.0'
import requests
import os
import subprocess as sp
import argparse
from typing import List


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
    
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    response = requests.get(url)
    if response.status_code != 200:
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
    
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Invalid URL: {url}")
    
    second_row = response.text.split("\n")[1]
    return second_row.split("\t")[1].split(";")

import signal, os

def handler(signum, frame):
    raise TimeoutError(f'Download timeout reached, trying again!')

def download_data(accession: str, urls: List[str]) -> None:
    """
    Download data from the ENA.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.
    urls : str
        The URLs of the data to download.

    Returns
    -------
    None
    """
    if not os.path.exists(accession):
        os.mkdir(accession)

    home = os.path.expanduser("~")
    ascp = os.path.join(home, '.aspera/cli/bin/ascp')
    opensshfile = os.path.join(home, '.aspera/cli/etc/asperaweb_id_dsa.openssh')

    i=0
    while i<3:
        signal.signal(signal.SIGALRM, handler)
        signal.alarm(2)
        try:
            for url in urls:
                path = url.replace('ftp.sra.ebi.ac.uk/', 'era-fasp@fasp.sra.ebi.ac.uk:')
                sp.run([
                    ascp, '-T', '-l', '300m', '-P', '33001', '-i', opensshfile, 
                    path, accession + '/'
                ], check=True)
        except:
            i+=1
            continue
    if i==3:
        raise TimeoutError(f"Download failed after 3 attempts")
    return None

def main(accession: str) -> None:
    """
    Function that calls all the other functions to download data from the ENA.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.

    Returns
    -------
    None
    """
    
    is_valid_accession(accession)
    
    paths = extract_data_path(accession)

    download_data(accession, paths)

    
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

    args = argparser.parse_args()

    main(args.accession)
