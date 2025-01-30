.. _sec-pre_imputation_qc:
.. _preimp_qc:

===================================
Pre-Imputation Quality Control (QC)
===================================

Detecting and correcting issues such as genotyping errors, sample handling errors, population stratification etc
is important in GWAS. The :code:`preimp_qc` module addresses these issues and cleans (QC) your data. Below is a flow diagram
of the filters applied when QC'ing input data:

.. image:: images/qc_workflow.png
   :width: 1000px
   :height: 1900px
   :scale: 50 %
   :align: center


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
   * - :code:`--liftover`
     - Liftover input data to GRCh38 or not, default is False. Running :code:`preimp_qc` with :code:`--liftover` will activate liftover
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


Examples
########

All the code below assumes the user already has a Dataproc cluster running as described in the `previous section <qb.html>`_

You can run pre-imputation qc using the :code:`preimp_qc` module (1) inside a python script; or (2) via the command line

1. Python script - submitting a python script to a cluster from local machine (Highly recommended)

- First create a python script on your local machine as below

    .. code-block:: python

        import gwaspy.preimp_qc as qc
        qc.preimp_qc.preimp_qc(dirname="gs://my-gcs/bucket/test_data/", basename="my_data_basename",
                               input_type="my_input_type")

- Then run the following command to submit the script to the Dataproc cluster named `my-cluster-name`

    .. code-block:: sh

        hailctl dataproc submit my-cluster-name qc_script.py

2. Command line - requires user to SSH'ed to a cluster

Users may encounter `this error <https://hail.zulipchat.com/#narrow/channel/128581-Cloud-support/topic/Running.20GWASpy.20on.20hailctl.20cluster.20-.20file.20not.20found.20exception>`_ when trying to run things from the command line

- This requires the user to be inside (`gcloud compute ssh`) the Dataproc cluster with GWASpy already installed

    .. code-block:: sh

        gcloud compute ssh "my-cluster-name-m"
        preimp_qc --dirname gs://my-gcs/bucket/test_data/ --basename my_data_basename --input-type my_input_type
