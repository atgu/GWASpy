#!/bin/bash

set -ex

PLATFORM="${OSTYPE}"

# install pylatex depedencies
case "$PLATFORM" in
    darwin*)
        install-pylatex-dependencies() {
            brew install --cask mactex
            eval "$(/usr/libexec/path_helper)"
        }
        install-git-lfs() {
            brew install git-lfs
            brew upgrade git-lfs
            git config --global lfs.batch false
            git lfs install
        }
        ;;
    linux*)
        install-pylatex-dependencies() {
            yes Y | apt-get install texlive-pictures texlive-science texlive-latex-extra latexmk
        }
        install-git-lfs() {
            curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
            yes Y | apt-get install git-lfs
        }
        ;;
    *)
        echo "unsupported platform $PLATFORM."
        ;;
esac

install-pylatex-dependencies
install-git-lfs
