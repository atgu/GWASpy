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

For usage, please visit [GWASpy](https://gwaspy.readthedocs.io/en/latest/index.html)