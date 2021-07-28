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
git clone https://github.com/atgu/GWASpy.git
cd GWASpy/
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
``--input-type`` | Input type, Hail MT, PLINK or VCF
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

**Argument** | **Description**
--- | ---
``--ref-dirname`` | Path to where reference data is
``--ref-basename`` | Reference basename 
``--ref-info`` | Path to reference information. Tab-delimited file with sample IDs and their SuperPop labels
``--reference`` | Genome reference build. Default is GRCh38
``--with-ref`` | Run PCA with or without reference. Default is without
``--data-dirname`` | Path to where the data is
``--data-basename`` | Data basename
``--input-type`` | Data input type. Choices are Hail MT, PLINK, and VCF
``--maf`` | include only SNPs with MAF >= NUM in PCA. Default is 0.05
``--hwe`` | include only SNPs with HWE >= NUM in PCA. Default is 1e-03 
``--geno`` | include only SNPs with call-rate > NUM. Default is 0.98
``--ld-cor`` | Squared correlation threshold (exclusive upper bound). Must be in the range [0.0, 1.0]. Default is 0.2
``--ld-window`` | Window size in base pairs (inclusive upper bound). Default is 250000
``--relatedness-method`` | Method to use for the inference of relatedness. Default is pc_relate
``--relatedness-thresh`` | Threshold value to use in relatedness checks. Default is 0.98
``--prob`` | Minimum probability of belonging to a given population for the population to be set. Default is 0.8
``--out-dir`` | Path to where output files will be saved

(3) Haplotype Phasing
---------------------
You can run haplotype phasing using Eagle_v2.4.1. We will ad support for SHAPEIT at a later stage

**Before running phasing, make sure that chromosomes Y and MT have been removed from the VCF, and ONLY chromosomes 1-22, X are present.**

**3.1. Running phasing without a reference panel:**
```bash
phasing --input-vcfs gs://path/to/vcf_files.txt --out-dir gs://path/to/output/directory
```

**3.2. Running phasing with a reference panel:**
```bash
phasing --input-vcfs gs://path/to/vcf_files.txt --vcf-ref gs://path/to/reference_panel.vcf --out-dir gs://path/to/output/directory
```
In the commands above, ``vcf_files.txt`` is a text file, with no header, containing path to each VCF. One VCF path per line

**Argument** | **Description**
--- | ---
``--input-vcfs`` | Path to where text file containing VCF(s) for target genotypes paths is
``--vcf-ref`` | VCF file for reference haplotypes
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
