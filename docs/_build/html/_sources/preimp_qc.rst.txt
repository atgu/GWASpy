.. _sec-pre_imputation_qc:

====================================
Pre-Imputation Quality Control (QC)
====================================

Detecting and correcting issues such as genotyping errors, sample handling errors, population stratification etc
is important in GWAS. The :code:`preimp_qc` module addresses these issues and cleans (QC) your data. Below is a flow diagram
of the filters applied when QC'ing input data:

.. image:: images/qc_workflow.png
   :width: 1000px
   :height: 1900px
   :scale: 50 %
   :align: center

Examples
########
You can run pre-imputation qc using the :code:`preimp_qc` module: (1) via the command line; or (2) inside a python script

#. Command line

    .. code-block:: sh

        preimp_qc --dirname data/ --basename sim_sim2a_eur_sa_merge.miss --input-type plink

#. Inside a python script

    .. code-block:: python

        import gwaspy.preimp_qc as qc
        qc.preimp_qc.preimp_qc(input_type="plink", dirname="data/", basename="sim_sim2a_eur_sa_merge.miss")


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
     - Input type. Options: [:code:`hail`, :code:`plink`, :code:`vcf`]
   * - :code:`--export-type`
     - Export type. Options: [:code:`hail`, :code:`plink`, :code:`vcf`]
   * - :code:`--out-dir`
     - Directory path to where output files are going to be saved
   * - :code:`--annotations`
     - Annotations file to be used for annotating sample with information such as Sex and Phenotype
   * - :code:`--reference`
     - Reference genome build. Default is GRCh38. Options: [:code:`GRCh37`, :code:`GRCh38`]
   * - :code:`--report`
     - Generate a QC PDF report or not. Default is True
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

Output(s)
##########
* QC'ed file(s) i.e. file with all the variants and/or samples that fail QC filters removed
* A detailed PDF QC report including pre- and post-QC variant/sample counts, figures such as Manhattan and QQ plots etc.
