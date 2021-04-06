preimp_qc
=============

A Python package for performing GWAS QC, PCA, and genotype imputation.

Installation
------------

For now you can install preimp_qc and its dependencies using the command below. In the near future, it will be uploaded to pypi

.. code:: bash

   bash env-setup.sh

Usage
-----

.. code:: bash

   $ preimp_qc --dirname data/ --basename sim_sim2a_eur_sa_merge.miss --input-type plink
   # in the example above, inside the directory data/, there will be 3 PLINK file sim_sim2a_eur_sa_merge.*{bed,bim,fam}

Arguments and Options
---------------------

+------------------------+--------------------------------------------+
| **Argument/ Option**   | **Description**                            |
+========================+============================================+
| ``--dirname``          | Path to where the data is                  |
+------------------------+--------------------------------------------+
| ``--basename``         | Data basename                              |
+------------------------+--------------------------------------------+
| ``--input-type``       | Input type, plink or vcf                   |
+------------------------+--------------------------------------------+
| ``--annotations``      | Annotations file to be used for annotating |
|                        | the VCF file (ONLY for VCF input)          |
+------------------------+--------------------------------------------+
| ``--reference``        | Reference genome build e.g. GRCh37, GRCh38 |
+------------------------+--------------------------------------------+
| ``--pre-geno``         | include only SNPs with missing-rate < NUM  |
|                        | (before ID filter), important for post     |
|                        | merge of multiple platforms                |
+------------------------+--------------------------------------------+
| ``--mind``             | include only IDs with missing-rate < NUM   |
+------------------------+--------------------------------------------+
| ``--fhet-aut``         | include only IDs within NUM < FHET < NUM   |
+------------------------+--------------------------------------------+
| ``--fhet-y``           | include only female IDs with fhet < NUM    |
+------------------------+--------------------------------------------+
| ``--fhet-x``           | include only male IDs with fhet > NUM      |
+------------------------+--------------------------------------------+
| ``--geno``             | include only SNPs with missing-rate < NUM  |
+------------------------+--------------------------------------------+
| ``--midi``             | include only SNPs with missing-rate        |
|                        | -difference ("case/control) < NUM          |
+------------------------+--------------------------------------------+
| ``--withpna``          | include monomorphic (invariant) SNPs       |
+------------------------+--------------------------------------------+
| ``--maf``              | include only SNPs with MAF >= NUM          |
+------------------------+--------------------------------------------+
| ``--hwe-th-con``       | HWE_controls < NUM                         |
+------------------------+--------------------------------------------+
| ``--hwe-th-cas``       | HWE_cases < NUM                            |
+------------------------+--------------------------------------------+