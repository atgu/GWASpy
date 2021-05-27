# GWASpy

<!-- badges: start -->
![install status](https://github.com/atgu/GWASpy/actions/workflows/install-ci.yml/badge.svg)

Genome-wide association studies pypeline (GWASpy): A Python package for performing GWAS QC, PCA, and genotype imputation.

## Installation


For now you can install GWASpy and its dependencies using the command below. In the near future, it will be uploaded to pypi

1. Install non-Python dependencies. The above command assumes the user already has [Homebrew](https://brew.sh/)
(MacOS) or [apt-get](https://linux.die.net/man/8/apt-get) (Linux) installed.

```bash
curl https://raw.githubusercontent.com/atgu/preimp_qc/main/env-setup.sh | bash
```

2. Clone the repo and install GWASpy and the required packages.
Once the package is on [pypi](https://pypi.org/), it won't be ncessary to run ``pip install -r requirements.txt``

```bash
# Install preimp_qc and its dependencies
git clone https://github.com/atgu/preimp_qc
cd preimp_qc/
pip3 install -r requirements.txt
python3 setup.py sdist
pip3 install dist/GWASpy-0.1.0.tar.gz
```

## Usage

## (1) Pre-imputation QC


You can run pre-imputation using the ``preimp_qc`` module (1) via the command line or (2) inside a python script

```bash
# (1) command line
$ preimp_qc --dirname data/ --basename sim_sim2a_eur_sa_merge.miss --input-type plink

# (2) inside a python script
import gwaspy.preimp_qc as qc
qc.preimp_qc.preimp_qc(input_type="plink", dirname="data/", basename="sim_sim2a_eur_sa_merge.miss")

# in the examples above, inside the directory data/, there will be 3 PLINK file sim_sim2a_eur_sa_merge.*{bed,bim,fam}
```

## Arguments and Options

**Argument** | **Description**
--- | ---
``--dirname`` | Path to where the data is
``--basename`` | Data basename
``--input-type`` | Input type, plink or vcf
``--annotations`` | Annotations file to be used for annotating<br>the VCF file (ONLY for VCF input)
``--reference`` | Reference genome build e.g. GRCh37, GRCh38
``--pre-geno`` | include only SNPs with missing-rate < NUM (before ID filter), important for post merge of multiple platforms
``--mind`` | include only IDs with missing-rate < NUM
``--fhet-aut`` | include only IDs within NUM < FHET < NUM
``--fhet-y`` | include only female IDs with fhet < NUM
``--fhet-x`` | include only male IDs with fhet > NUM
``--geno`` | include only SNPs with missing-rate < NUM
``--midi`` | include only SNPs with missing-rate-difference ("case/control) < NUM
``--withpna`` | include monomorphic (invariant) SNPs
``--maf`` | include only SNPs with MAF >= NUM
``--hwe-th-con`` | HWE_controls < NUM
``--hwe-th-cas`` | HWE_cases < NUM

(2) Principal Component Analysis
--------------------------------
GWASpy has two way to run PCA: (1) without a reference; or (2) with a reference (Coming). The output result is a tsv file with scores for each sample and a PDF file with PCA plots

Below is an example of how to run PCA without a reference using the example data ``data/1kg_annotated.mt``

```bash
pca --dirname data/ --basename 1kg_annotated --out-dir data/ --input-type hail --reference grch37
```

(3) Haplotype Phasing
---------------------
You can run haplotype phasing using Eagle_v2.4.1. We will ad support for SHAPEIT at a later stage

Here's how you can run phasing:
```bash
phasing --input-vcfs gs://path/to/vcf_files.txt --out-dir gs://path/to/output/directory
```
In the command above, ``vcf_files.txt`` is a text file, with no header, containing path to each VCF. One VCF path per line

**Argument** | **Description**
--- | ---
``--input-vcfs`` | Path to where text file containing VCF(s) paths is
``--local`` | Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud 
``--software`` | Software to use for phasing. Default is Eagle
``--cpu`` | Number of CPUs to use. Default is 8
``--memory`` | Memory to use. Default is standard which correspond to ~4Gi/core. lowmem ~1Gi/core and highmem ~7Gi/core
``--storage`` | Storage to use for the job in gigabytes. Default is 50 Gi
``--threads`` | Number of threads to use in phasing. Default is 16
``--out-dir`` | Path to where output files will be saved

(4) Genotype imputation
--------------------------------
Coming...
