*Tracking COVID-19 Variants*
---------------------------------------------------------------

Video tutorial
-------------
If you have a WHO account, you should be able to
access
[this video tutorial](https://worldhealthorg-my.sharepoint.com/:v:/g/personal/campbellf_who_int/EZu9wkwnGbhLl_iIsTMDD8gBe_j9At-qkm7EkqQ22yvrVQ?e=C6wg1h&isSPOFile=1),
which will walk you through the analysis. If
this doesn't work, contact me under campbellf@who.int or follow the brief guide
below.


Downloading data
-------------

- Download the "Genomic epidemiology" `metadata` database from
the [GISAID](https://www.epicov.org/epi3/frontend#278a01) website (under the
"Downloads" tab). You will need a GISAID account for this.
- Extract the `metadata.tsv` file and save it under the `\data` folder in the
  format `genomic_epi_YYYY-MM-DD.tsv`. 

Running the analyis
-------------

- Open the `covariant_tracking.RProj` file if using RStudio, otherwise make sure
  your working directory is somewhere within this folder
- Open the `run_analysis.R` file, which will guide you through the analysis.
