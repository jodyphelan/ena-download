# ENA-Download

You can use the this package to download data from ENA. The CLI is installed with the package. To use it, you need to have a aspera client installed. To find out more about the aspera client, visit [this page](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html#using-aspera). Make sure you install it in your home holder at `~/.aspera`.

To download just provide the accession number of the data you want to download. A folder will be created with the accession number as the name and the data will be downloaded to that folder.

```
ena_download <acession>
```

