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
        ;;
    linux*)
        install-pylatex-dependencies() {
            yes Y | apt-get install texlive-pictures texlive-science texlive-latex-extra latexmk
        }
        ;;
    *)
        echo "unsupported platform $PLATFORM."
        ;;
esac

install-pylatex-dependencies
