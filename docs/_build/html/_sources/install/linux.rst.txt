===========================
Install GWASpy on GNU/Linux
===========================

- Install Java 8.
- Install Python 3.6+.
- Install a recent version of the C and C++ standard libraries. GCC 5.0, LLVM
  version 3.4, or any later versions suffice.
- Install BLAS and LAPACK.
- Install TeX Live
- Install GWASpy using pip.

On a recent Debian-like system, the following should suffice:

.. code-block:: sh

   apt-get install -y \
       openjdk-8-jre-headless \
       g++ \
       python3.6 python3-pip \
       libopenblas-base liblapack3 \
       texlive-pictures texlive-science texlive-latex-extra latexmk
   python3.6 -m pip install gwaspy
