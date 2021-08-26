.. _sec-phasing:

===================
Haplotype Phasing
===================

Knowing the phase of a haplotype can allow us to impute low frequency variants, this makes haplotype phasing an
important step before genotype imputation. GWASpy has a module, :code:`phasing`, for performing phasing. Phasing can
be run with or without a reference panel using either Eagle2 or SHAPEIT4

.. note::
    If the data set you wish to phase contains more than twice as many samples as the largest reference panel
    available to you, then using a reference panel is unlikely to give much of a boost in phasing accuracy.

.. toctree::
   :maxdepth: 1

        Phasing without a reference panel <phasing/without_ref.rst>
        Phasing with a reference panel <phasing/with_ref.rst>

Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--input-vcfs`
     - Path to where text file containing VCF(s) for target genotypes paths is
   * - :code:`--vcf-ref`
     - VCF file for reference haplotypes if phasing with a reference panel
   * - :code:`--local`
     - Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud
   * - :code:`--billing-project`
     - Billing project to be used for the job(s)
   * - :code:`--bucket`
     - Bucket associated with the billing project
   * - :code:`--software`
     - Software to use for phasing. Options: [:code:`eagle`, :code:`shapeit`]. Default is Eagle
   * - :code:`--cpu`
     - Number of CPUs to use. Default is 8
   * - :code:`--memory`
     - Memory to use. Default is :code:`standard` which correspond to ~4Gi/core. :code:`lowmem` ~1Gi/core and :code:`highmem` ~7Gi/core
   * - :code:`--storage`
     - Storage to use for the job in gigabytes. Default is 50 Gi
   * - :code:`--threads`
     - Number of threads to use in phasing. Default is 16
   * - :code:`--out-dir`
     - Path to where output files will be saved

Output
######
For both Eagle and SHAPEIT, the resulting output is a VCF file per chromosome with phased haplotypes.


.. note::
    By default, Eagle will output a VCF file with phased GT and other fields that were in the unphased VCF, whereas
    SHAPEIT will ONLY output the GT field. This will result in phased files generated using Eagle being bigger in size
    than those generate using SHAPEIT.
