.. _sec-pre_imputation_qc:

===================================
Pre-Imputation Quality Control (QC)
===================================

Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--dirname`
     - Path to where the data is
   * - :code:`--basename`
     - Data basename
   * - :code:`--input-type`
     - Input type, Hail MT, PLINK or VCF
   * - :code:`--annotations`
     - Annotations file to be used for annotating<br>the VCF file (ONLY for VCF input)
   * - :code:`--reference`
     - Reference genome build e.g. GRCh37, GRCh38
   * - :code:`--pre-geno`
     - include only SNPs with missing-rate < NUM (before ID filter), important for post merge of multiple platforms
   * - :code:`--mind`
     - include only IDs with missing-rate < NUM
   * - :code:`--fhet-aut`
     - include only IDs within NUM < FHET < NUM
   * - :code:`--fstat-y`
     - include only female IDs with fhet < NUM
   * - :code:`--fstat-x`
     - include only male IDs with fhet > NUM
   * - :code:`--geno`
     - include only SNPs with missing-rate < NUM
   * - :code:`--midi`
     - include only SNPs with missing-rate-difference (case/control) < NUM
   * - :code:`--withpna`
     - include monomorphic (invariant) SNPs
   * - :code:`--maf`
     - include only SNPs with MAF >= NUM
   * - :code:`--hwe-th-con`
     - HWE_controls < NUM
   * - :code:`--hwe-th-cas`
     - HWE_cases < NUM
