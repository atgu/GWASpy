#!/bin/bash

set -ex

# Install preimp_qc and its dependencies
git clone https://github.com/atgu/preimp_qc
cd preimp_qc/
pip install -r requirements.txt
python setup.py sdist
pip install dist/preimp_qc-0.1.0.tar.gz

PLATFORM="${OSTYPE}"

# install pylatex depedencies
case "$PLATFORM" in
    darwin*)
        install-pylatex-dependencies() {
            brew install --cask mactex
            eval "$(/usr/libexec/path_helper)"
        }
        ;;
    linux*)
        install-pylatex-dependencies() {
            sudo apt-get install texlive-pictures texlive-science texlive-latex-extra latexmk
        }
        ;;
    *)
        echo "unsupported platform $PLATFORM."
        ;;
esac

install-pylatex-dependencies
