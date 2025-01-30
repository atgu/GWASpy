.. _sec-installation:

=================
Installing GWASpy
=================

GWASpy leverages Hail to enable efficient processing of data directly from Google Cloud. As such, the first step is to
install Hail as per instructions `here <https://hail.is/docs/0.2/install/macosx.html>`_. After you have installed Hail, GWASpy can be easily installed using

.. code-block:: sh

   pip install gwaspy

It is important to note that the command above will install GWASpy locally (or wherever you ran the command). For the
:code:`phasing` and :code:`imputation` modules using Hail Batch, this is enough. For the :code:`preimp_qc` and
:code:`pca` modules using Hail Query, however, you have to ensure that the dataproc cluster has GWASpy, and there are
examples showing how to do this in the :ref:`preimp_qc` and :ref:`pca` sections.
