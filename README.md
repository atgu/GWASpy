# GWASpy

<!-- badges: start -->
![install status](https://github.com/atgu/GWASpy/actions/workflows/install-ci.yml/badge.svg)
![pca status](https://github.com/atgu/GWASpy/actions/workflows/pca-ci.yml/badge.svg)

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
pip install -r requirements.txt
python setup.py sdist
pip install dist/GWASpy-0.1.0.tar.gz
```

## Usage

## (1) preimp_qc


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
Coming...

(3) Genotype imputation
--------------------------------
Coming...
