"""
This package facilitates the download of data from the ENA in fastq format.
To use it, you need to provide the accession number of the data you want to download.
"""
__version__ = '0.3.0'
import requests
import os
import subprocess as sp
import argparse
from typing import List, Dict
import sys
import json
import logging
from ftplib import FTP
import tempfile 
import shutil
import hashlib


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

def extract_data_path(accession: str) -> Dict[str, str]:
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
    
    md5_list = []
    for d in data:
        md5_list += d['fastq_md5'].split(";")
    
    md5s = dict(zip(files, md5_list))
    
    return md5s

import signal, os

def handler(signum, frame):
    """Signal handler for the download timeout."""
    raise TimeoutError(f'Download timeout reached, trying again!')

def md5sum(file: str) -> str:
    """
    Calculate the md5 checksum of a file.

    Parameters
    ----------
    file : str
        The path to the file to calculate the md5 checksum for.

    Returns
    -------
    str
        The md5 checksum of the file.
    """

    hash_md5 = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def ftp_download_data(accession: str, output_directory: str, files: Dict[str, str]) -> None:
    """
    Download data from the ENA.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.
    output_directory : str
        The directory to download the data to.
    urls : str
        The URLs of the data to download.
    md5s : dict
        The md5 checksums of the files to download.

    Returns
    -------
    None
    """

    ftp = FTP('ftp.sra.ebi.ac.uk')
    ftp.login('anonymous')

    urls = list(files.keys())
    md5s = files

    logging.debug(f"Downloading data for {accession}")

    with tempfile.TemporaryDirectory() as tmpdirname:
        for url in urls:
            logging.debug(f"Downloading {url} into {tmpdirname}")
            location = url.replace('ftp.sra.ebi.ac.uk', '')
            filename = url.split('/')[-1]
            ftp.retrbinary(f'RETR {location}', open(os.path.join(tmpdirname, filename), 'wb').write)

            # Check md5 checksum
            md5 = md5s[url]
            download_md5 = md5sum(os.path.join(tmpdirname, filename))
            md5_match = md5 == download_md5
            logging.debug(f"MD5 checksum for {filename} = {md5_match}")
            if not md5_match:
                raise ValueError(f"MD5 checksum failed for {url}")

        if accession.startswith('SAM'):
            forward_reads = sorted([f'{tmpdirname}/{f}' for f in os.listdir(tmpdirname) if f.endswith('_1.fastq.gz')])
            reverse_reads = sorted([f'{tmpdirname}/{f}' for f in os.listdir(tmpdirname) if f.endswith('_2.fastq.gz')])
            if len(forward_reads) == 0 or len(reverse_reads) == 0:
                raise ValueError(f"Download failed for {accession}")

            sp.run(f"cat {' '.join(forward_reads)} > {os.path.join(tmpdirname, accession + '_1.fastq.gz')}", shell=True, check=True)
            sp.run(f"cat {' '.join(reverse_reads)} > {os.path.join(tmpdirname, accession + '_2.fastq.gz')}", shell=True, check=True)

            # remove the original files
            for f in forward_reads + reverse_reads:
                os.remove(os.path.join(tmpdirname, f))
        
        # move the files to the output directory
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)
        for f in os.listdir(tmpdirname):
            logging.debug(f"Moving {f} to {output_directory}")
            shutil.move(os.path.join(tmpdirname, f), os.path.join(output_directory, f))

    return None

def main(accession: str,  output_directory: str) -> None:
    """
    Function that calls all the other functions to download data from the ENA.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.
    mode : str
        The mode of download: ftp or ascp.
    timeout : int
        The timeout in seconds for the download to complete. Default is 300 seconds.

    Returns
    -------
    None
    """
    
    is_valid_accession(accession)
    
    files = extract_data_path(accession)


    ftp_download_data(
        accession=accession,
        output_directory=output_directory,
        files=files
    )

    
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
    argparser.add_argument('--outdir', default=".", type=str, help='Output directory to download the data to')
    argparser.add_argument('--debug', action='store_true', help='Print debug information')
    args = argparser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    main(args.accession,args.outdir)
