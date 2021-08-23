.. _sec-phasing:

===================================
Haplotype Phasing
===================================

Write something about haplotype phasing and include examples

.. toctree::
   :maxdepth: 1

        Phasing without a reference panel <phasing/without_ref.rst>
        Phasing with a reference panel <phasing/with_ref.rst>

Arguments and options
#####################

.. list-table::
   :widths: 25 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--input-vcfs`
     - Path to where text file containing VCF(s) for target genotypes paths is
   * - :code:`--vcf-ref`
     - VCF file for reference haplotypes
   * - :code:`--local`
     - Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud
   * - :code:`--software`
     - Software to use for phasing. Default is Eagle
   * - :code:`--cpu`
     - Number of CPUs to use. Default is 8
   * - :code:`--memory`
     - Memory to use. Default is standard which correspond to ~4Gi/core. lowmem ~1Gi/core and highmem ~7Gi/core
   * - :code:`--storage`
     - Storage to use for the job in gigabytes. Default is 50 Gi
   * - :code:`--threads`
     - Number of threads to use in phasing. Default is 16
   * - :code:`--out-dir`
     - Path to where output files will be saved
