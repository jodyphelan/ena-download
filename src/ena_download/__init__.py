"""
This package facilitates the download of data from the ENA in fastq format.
To use it, you need to provide the accession number of the data you want to download.
"""
__version__ = '0.3.1'
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
from tqdm import tqdm
import time


logging.basicConfig(level=logging.INFO)

# Global flag to control progress bar display
_show_progress = True

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
    time.sleep(1)
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    parameters = {
        "result": "read_run",
        "includeAccessions": accession,
        "limit": "10",
        "format": "json"
    }
    response = requests.get(url, params=parameters)
    
    data = json.loads(response.text)

    if len(data) == 0:
        raise ValueError(f"Invalid accession number: {accession}")
    return True

def get_accession_info(accession: str) -> Dict[str, str]:
    """
    Get information about the accession number.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.

    Returns
    -------
    dict
        A dictionary containing information about the accession number.

    Examples
    --------
    >>> get_accession_info("ERR11466368")
    {'run_accession': 'ERR11466368', 'instrument_model': 'Illumina NovaSeq 6000', 'instrument_platform': 'ILLUMINA', 'library_strategy': 'WGS', 'library_layout': 'PAIRED', 'library_source': 'GENOMIC'}
    """
    time.sleep(1)
    logging.debug(f"Getting accession info for {accession}")
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    parameters = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession,instrument_model,instrument_platform,library_strategy,library_layout,library_source",
        "format": "json"
    }

    response = requests.get(url, params=parameters)
    if response.status_code != 200:
        raise ValueError(f"Invalid URL: {url}")
    
    logging.debug(f"Response: {response.text}")
    data = json.loads(response.text)
    logging.debug(f"Data found for {accession}: {data}")

    if len(data) == 0:
        raise ValueError(f"No data found for {accession}")

    return data

def extract_data_path(accession: str, platform: str, library_strategy: str) -> Dict[str, str]:
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
    data = get_accession_info(accession)
    files = {}
    for d in data:
        if d['instrument_platform'] == platform and d['library_strategy'] == library_strategy:
            tmpfiles = d['fastq_ftp'].split(";")
            tmpmd5s = d['fastq_md5'].split(";")
            for f, m in zip(tmpfiles, tmpmd5s):
                files[f] = m
        else:
            logging.debug(f"Skipping {d['run_accession']} due to platform/library strategy mismatch (found: {d['instrument_platform']}/{d['library_strategy']}, expected: {platform}/{library_strategy})")
    
    if len(files) == 0:
        raise ValueError(f"No data found for {accession}")
    
    return files

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

def http_get_file(url: str, tmpdirname: str) -> None:
    """Download a single file from ENA via HTTP with up to 3 retries.

    Parameters
    ----------
    url : str
        FTP-style URL to the file (e.g., ftp.sra.ebi.ac.uk/vol1/fastq/.../file.fastq.gz)
        Will be converted to HTTPS URL automatically.
    tmpdirname : str
        Temporary directory path to write the file during download
    """
    max_attempts = 3
    backoff_base = 1  # seconds

    # Convert FTP URL to HTTPS
    if url.startswith('ftp.sra.ebi.ac.uk'):
        http_url = f"https://{url}"
    else:
        http_url = url

    filename = url.split('/')[-1]
    dest_path = os.path.join(tmpdirname, filename)
    logging.debug(f"Downloading {http_url} into {tmpdirname} via HTTP")

    for attempt in range(1, max_attempts + 1):
        time.sleep(1)  # small wait before each attempt
        # Ensure any partial file from a previous attempt doesn't remain
        if os.path.exists(dest_path):
            try:
                os.remove(dest_path)
            except Exception:
                pass

        try:
            response = requests.get(http_url, stream=True, timeout=300)
            response.raise_for_status()

            # Get file size if available
            total_size = int(response.headers.get('content-length', 0))

            with open(dest_path, 'wb') as f:
                if total_size and _show_progress:
                    with tqdm(total=total_size, unit='B', unit_scale=True, desc=filename) as pbar:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
                else:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)

            # If we reached here without raising, download succeeded
            return
        except Exception as e:
            # Log and retry with exponential backoff
            if attempt < max_attempts:
                wait = backoff_base * (2 ** (attempt - 1))
                logging.warning(f"Attempt {attempt} to download {filename} via HTTP failed: {e}. Retrying in {wait}s...")
                time.sleep(wait)
            else:
                logging.error(f"Failed to download {filename} via HTTP after {max_attempts} attempts.")
                raise

def ftp_get_file(ftp: FTP, url: str, tmpdirname: str) -> None:
    """Download a single file from ENA FTP with up to 3 retries.

    Parameters
    ----------
    ftp : FTP
        An active ftplib.FTP connection to ftp.sra.ebi.ac.uk
    url : str
        Full ftp URL to the file (e.g., ftp.sra.ebi.ac.uk/vol1/fastq/.../file.fastq.gz)
    tmpdirname : str
        Temporary directory path to write the file during download
    """
    max_attempts = 3
    backoff_base = 1  # seconds


    logging.debug(f"Downloading {url} into {tmpdirname}")
    location = url.replace('ftp.sra.ebi.ac.uk', '')
    filename = url.split('/')[-1]

    for attempt in range(1, max_attempts + 1):
        time.sleep(1)  # small wait before each attempt
        # Ensure any partial file from a previous attempt doesn't remain
        dest_path = os.path.join(tmpdirname, filename)
        if os.path.exists(dest_path):
            try:
                os.remove(dest_path)
            except Exception:
                # If removal fails, we'll still overwrite in 'wb' mode
                pass

        try:
            # Try to get the size for progress reporting
            try:
                total_size = ftp.size(location)
            except Exception:
                logging.info(f"Error getting size for {location}... not reporting progress.")
                total_size = None

            with open(dest_path, 'wb') as f:
                if total_size and _show_progress:
                    with tqdm(total=total_size, unit='B', unit_scale=True, desc=filename) as pbar:
                        def callback(data):
                            f.write(data)
                            pbar.update(len(data))
                        ftp.retrbinary(f'RETR {location}', callback)
                else:
                    # Fallback without progress bar
                    def write_chunk(data):
                        f.write(data)
                    ftp.retrbinary(f'RETR {location}', write_chunk)

            # If we reached here without raising, download succeeded
            return
        except Exception as e:
            # Log and retry with exponential backoff
            if attempt < max_attempts:
                wait = backoff_base * (2 ** (attempt - 1))
                logging.warning(f"Attempt {attempt} to download {filename} failed: {e}. Retrying in {wait}s...")
                time.sleep(wait)
                # best-effort: keep the FTP connection alive between retries
                try:
                    ftp.voidcmd('NOOP')
                except Exception:
                    # If NOOP fails, continue anyway; caller opened the connection
                    pass
            else:
                logging.error(f"Failed to download {filename} after {max_attempts} attempts.")
                raise

def ftp_download_data(accession: str, output_directory: str, files: Dict[str, str]) -> None:
    """Download data from the ENA via FTP, with HTTP fallback.

    Parameters
    ----------
    accession : str
        The accession number of the data to download.
    output_directory : str
        The directory to download the data to.
    files : dict
        Dict mapping URLs to their md5 checksums.

    Returns
    -------
    None
    """

    urls = list(files.keys())
    md5s = files

    logging.debug(f"Downloading data for {accession}")

    with tempfile.TemporaryDirectory() as tmpdirname:
        for url in urls:
            filename = url.split('/')[-1]
            download_success = False

            # Try FTP first
            try:
                ftp = FTP('ftp.sra.ebi.ac.uk')
                ftp.login('anonymous')
                ftp_get_file(ftp, url, tmpdirname)
                ftp.quit()
                download_success = True
                logging.info(f"Successfully downloaded {filename} via FTP")
            except Exception as ftp_error:
                logging.warning(f"FTP download failed for {filename}: {ftp_error}")
                try:
                    ftp.quit()
                except:
                    pass

                # Fallback to HTTP
                logging.info(f"Attempting HTTP download for {filename}...")
                try:
                    http_get_file(url, tmpdirname)
                    download_success = True
                    logging.info(f"Successfully downloaded {filename} via HTTP")
                except Exception as http_error:
                    logging.error(f"HTTP download also failed for {filename}: {http_error}")
                    raise ValueError(f"Failed to download {filename} via both FTP and HTTP") from http_error

            if not download_success:
                raise ValueError(f"Failed to download {filename}")

            # Check md5 checksum
            md5 = md5s[url]
            download_md5 = md5sum(os.path.join(tmpdirname, filename))
            md5_match = md5 == download_md5
            logging.debug(f"MD5 checksum for {filename} = {md5_match} (expected: {md5}, got: {download_md5})")
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

def main(
        accession: str,  
        output_directory: str,
        platform: str,
        library_strategy: str,
        show_progress: bool = True
    ) -> None:
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
    show_progress : bool
        Whether to display progress bars during download. Default is True.

    Returns
    -------
    None
    """
    global _show_progress
    _show_progress = show_progress
    
    is_valid_accession(accession)
    
    files = extract_data_path(accession, platform, library_strategy)


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
    argparser.add_argument('--platform', type=str, default="ILLUMINA", help='Instrument platform to filter the data')
    argparser.add_argument('--library_strategy', type=str, default='WGS', help='Library strategy to filter the data')
    argparser.add_argument('--no-progress', action='store_true', help='Disable progress bar display')
    argparser.add_argument('--debug', action='store_true', help='Print debug information')
    argparser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = argparser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    main(args.accession, args.outdir, args.platform, args.library_strategy, show_progress=not args.no_progress)
