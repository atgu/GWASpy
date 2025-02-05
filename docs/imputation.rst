.. _sec-imputation:

===================
Genotype Imputation
===================

Genotype imputation is a process of estimating missing genotypes from the haplotype or genotype reference panel. It
allows you to accurately evaluate the evidence for association at genetic markers that are not directly genotyped.
GWASpy has a module, :code:`imputation`, for running imputation using IMPUTE5. Because imputation can be a computationally
intensive task, we run it on multiple chunks in parallel, then merge the imputed chunks together at the end. Below are
examples of how to run imputation using either the HGDP+1kGP or your own reference panel.

Examples
########

**1. HGDP+1kGP reference panel**

    .. code-block:: sh

        imputation --input-file gs://path/to/file.vcf.bgz --vcf-ref hgdp1kgp --output-filename my_outfilename --out-dir gs://path/to/output/dir --n-samples 1989 --n-ref-samples 4091 --billing-project my-billing-project

**2. Own reference panel**

    .. code-block:: python

        imputation --input-file gs://path/to/file.vcf.bgz --vcf-ref gs://path/to/ref_panel/ALL.chrCNUMBER.vcf --output-filename my_outfilename --out-dir gs://path/to/output/dir --n-samples 1989 --n-ref-samples 4091 --billing-project my-billing-project

.. warning::
    When using your own reference panel, make sure that you use the CNUMBER placeholder in the filename passed to --vcf-ref

Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--input-file`
     - Path to where the VCF or TSV with target VCF/BAM files is
   * - :code:`--vcf-ref`
     - Reference panel file to use for imputation
   * - :code:`--chromosomes`
     - Chromosome(s) to run imputation for. Default is :code:`all`
   * - :code:`--local`
     - Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud
   * - :code:`--billing-project`
     - Billing project to be used for the jobs
   * - :code:`--n-samples`
     - Number of target samples to be imputed. We use this to estimate resources for some of the jobs
   * - :code:`--n-ref-samples`
     - Number of reference samples. We use this to estimate resources for some of the jobs
   * - :code:`--software`
     - Software to use for phasing. Options: [:code:`beagle5`, :code:`impute5`]. Default is :code:`impute5`
   * - :code:`--output-filename`
     - Output filename without file extension
   * - :code:`--out-dir`
     - Path to where output files will be saved

Output
######
The resulting output is a VCF file per chromosome with imputed genotypes.
