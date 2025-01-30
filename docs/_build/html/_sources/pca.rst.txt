.. _sec-pca:
.. _pca:

============================
Principal Component Analysis
============================

Principal components analysis (PCA) can be used to detect and quantify the genetic structure of populations.
In GWASpy, the :code:`pca` module can be run in 3 different ways: (1) normal PCA without a reference panel; (2) joint PCA; or (3) Projection PCA.

.. toctree::
   :maxdepth: 1

        Normal PCA <pca/normal.rst>
        Joint PCA <pca/joint.rst>
        Projection PCA <pca/project.rst>

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
   * - :code:`--pca-type`
     - Type of PCA to run. Default is normal. Options: [:code:`normal`, :code:`project`, :code:`joint`]
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
   * - :code:`--npcs`
     - Number of PCs to use. Default is 20
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
