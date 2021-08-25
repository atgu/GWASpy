.. _sec-pca:

============================
Principal Component Analysis
============================

Principal components analysis (PCA) can be used to detect and quantify the genetic structure of populations.
In GWASpy, the :code:`pca` module can be run in two different ways: (1) without a reference panel; and (2) with a reference panel.

.. toctree::
   :maxdepth: 1

        PCA without a reference panel <pca/without_ref.rst>
        PCA with a reference panel <pca/with_ref.rst>

Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--ref-dirname`
     - Path to where reference data is
   * - :code:`--ref-basename`
     - Reference basename
   * - :code:`--ref-info`
     - Path to reference information. Tab-delimited file with sample IDs and their SuperPop labels
   * - :code:`--reference`
     - Genome reference build. Default is GRCh38. Options: [:code:`GRCh37`, :code:`GRCh38`]
   * - :code:`--with-ref`
     - Run PCA with or without reference. Default is without
   * - :code:`--data-dirname`
     - Path to where the data is
   * - :code:`--data-basename`
     - Data basename
   * - :code:`--input-type`
     - Data input type. Options: [:code:`hail`, :code:`plink`, :code:`vcf`]
   * - :code:`--maf`
     - include only SNPs with MAF >= NUM in PCA. Default is 0.05
   * - :code:`--hwe`
     - include only SNPs with HWE >= NUM in PCA. Default is 1e-03
   * - :code:`--geno`
     - include only SNPs with call-rate > NUM. Default is 0.98
   * - :code:`--ld-cor`
     - Squared correlation threshold (exclusive upper bound). Must be in the range [0.0, 1.0]. Default is 0.2
   * - :code:`--ld-window`
     - Window size in base pairs (inclusive upper bound). Default is 250000
   * - :code:`--relatedness-method`
     - Method to use for the inference of relatedness. Default is pc_relate. Options: [:code:`pc_relate`, :code:`ibd`, :code:`king`]
   * - :code:`--relatedness-thresh`
     - Threshold value to use in relatedness checks. Default is 0.98
   * - :code:`--prob`
     - Minimum probability of belonging to a given population for the population to be set. Default is 0.8
   * - :code:`--out-dir`
     - Path to where output files will be saved

Output
######
A tab-delimited file with the first 20 principal components (PCs)  computed and
graphical visualizations of the PCs are generated.
