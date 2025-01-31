.. _sec-tutorial:

========
Tutorial
========

This is a short tutorial on how to use the different modules of GWASpy.

1. Datasets
###########

We will be using simulated test data (on GRCh37) from RICOPILI. Below is how it can be downloaded and copied to a Google bucket

    .. code-block:: sh

        wget https://personal.broadinstitute.org/sawasthi/share_links/UzoZK7Yfd7nTzIxHamCh1rSOiIOSdj_gwas-qcerrors.py/sim_sim1a_eur_sa_merge.miss.{bed,bim,fam} .
        gsutil cp sim_sim1a_eur_sa_merge.miss.{bed,bim,fam} gs://my-gcs/bucket/test_data

2. Start a dataproc cluster with GWASpy installed
#################################################

The code below will start a cluster with GWASpy automatically installed. This fetches the GWASpy version on PyPI. You
can also install the GitHub version by replacing :code:`gwaspy` with :code:`git+https://github.com/atgu/GWASpy.git`

    .. code-block:: sh

            hailctl dataproc start gwaspy-tut --region=us-central1 --packages gwaspy --max-age 4h

3. Pre-imputation QC
####################

Next, we will QC the data using the default arguments in :code:`preimp_qc`. Since the reference panel we will be using
for phasing and imputation is on GRCh38, we also add a :code:`liftover` argument to liftover our input data to GRCh38. We
also set :code:`export_type` to :code:`vcf` because we will use the QC'ed file as input to phasing and imputation.

- First create a python script on your local machine as below

    .. code-block:: python

        import gwaspy.preimp_qc as qc
        qc.preimp_qc.preimp_qc(dirname="gs://my-gcs/bucket/test_data/", basename="sim_sim1a_eur_sa_merge.miss",
                               input_type="plink", reference="GRCh37", liftover=True, export_type="vcf")

- Then run the following command to submit the script to the Dataproc cluster named `gwaspy-tut`

    .. code-block:: sh

        hailctl dataproc submit gwaspy-tut qc_script.py

4. PCA
######

Using the QC'ed data from step 3 above, we will now run PCA.

- First create a python script on your local machine as below

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="gs://my-gcs/bucket/test_data/GWASpy/Preimp_QC",
                    data_basename="sim_sim1a_eur_sa_merge.miss_qced", out_dir="gs://my-gcs/bucket/test_data/",
                    input_type="vcf", reference="GRCh38", pca_type="normal")

- Then run the following command to submit the script to the Dataproc cluster named `gwaspy-tut`

    .. code-block:: sh

        hailctl dataproc submit gwaspy-tut pca_script.py

If you have real data, you can use (1) :code:`pca_type="project"` which will train a Random Forest model on the HGDP+1kGP reference and
use the trained model to classify samples in your data; or (2) :code:`pca_type="joint"` which will first find an intersection (variants) between
the HGDP+1kGP reference and your input, use the intersected HGDP+1kGP to train a RF model, then classify your input data. If
your data has a lot of variants (+million), :code:`pca_type="project"` usually gives plausible results. Otherwise you can try :code:`pca_type="joint"`

.. note::
    If you are a Broad user with Hail Batch access, you have to have python and GWASpy installed locally to be able to run
    phasing and imputation. For non-Broad users, we provide a nextflow implementation and the only thing you are required
    to do is have nextflow locally (nextflow executable file) and necessary permissions `as mentioned <qb.html>`_

5. Phasing and Imputation
#########################

**5.1 Hail Batch**

5.1.1. Phasing (should be ~$2 and takes ~40 minutes)

The example below is for running phasing, without a reference panel. If you want to use the HGDP+1kGP reference panel or
your own, simply add the :code:`--vcf-ref` argument `as explained here <phasing.html>`_

    .. code-block:: sh

        phasing --input-vcf gs://my-gcs/bucket/test_data/GWASpy/Preimp_QC/sim_sim1a_eur_sa_merge.miss_qced.vcf.bgz \
        --output-filename sim_sim1a_eur_sa_merge.miss_qced.phased --out-dir gs://my-gcs/bucket/test_data/GWASpy/phasing \
        --fill-tags --genome-build GRCh38 --billing-project my-billing-project

5.1.2. Imputation using IMPUTE5 (should be ~$4 and takes <20 minutes)

The example below is for running phasing, without a reference panel. If you want to use the HGDP+1kGP reference panel or
your own, simply add the :code:`--vcf-ref` argument `as explained here <phasing.html>`_

    .. code-block:: sh

        imputation --input-file gs://my-gcs/bucket/test_data/GWASpy/phasing/shapeit5/phase_common/sim_sim1a_eur_sa_merge.miss_qced.phased_chrCNUMBER.array.shapeit5_common.bcf \
        --vcf-ref hgdp1kgp --output-filename sim_sim1a_eur_sa_merge.miss_qced.phased.imputed --out-dir gs://my-gcs/bucket/test_data/GWASpy/imputation \
        --n-samples 1989 --n-ref-samples 4091 --billing-project my-billing-project

.. note::
    You may need to add :code:`HAIL_GENETICS_HAIL_IMAGE=hailgenetics/python-dill:3.9-slim` in front of the :code:`phasing`
    and :code:`imputation` commands if you are using a Python version other than 3.9, 3.10, or 3.11

**5.2. Nextflow**

Before we run the nextflow pipeline, you have to first download the following files and copy them to your bucket:
(1) common chunks and rare chunks files used to parallelize imputation across genomic regions; (2) genetic map files. SHAPEIT5 repo
has `chunks files <https://github.com/odelaneau/shapeit5/tree/main/resources/chunks/b38>`_ and `genetic map files <https://github.com/odelaneau/shapeit5/tree/main/resources/maps>`_.

Once you have the files on a Google bucket, you can update the :code:`params.json` file. Specifically, the things you need
to update are: :code:`input_vcf`, :code:`output_filename`, :code:`out_dir`, :code:`data_type`, :code:`common_chunks`,
:code:`rare_chunks`, :code:`genetic_maps`. If you have one input file per chromosome, set :code:`input_split_by_chrom` to :code:`true`

    .. code-block:: sh

        {
            "input_vcf": "gs://my-gcs/bucket/test_data/GWASpy/Preimp_QC/sim_sim1a_eur_sa_merge.miss_qced.vcf",
            "output_filename": "sim_sim1a_eur_sa_merge.miss_qced",
            "out_dir": "gs://my-gcs/bucket/test_data/GWASpy/nf_phase_impute",
            "impute": true,
            "fill_tags": true,
            "input_split_by_chrom": false,
            "vcf_ref": "gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_chrCNUMBER.filtered.SNV_INDEL.phased.shapeit5",
            "ref_format": "vcf",
            "data_type": "array", // or wgs
            "maf": 0.001,
            "common_chunks": "gs://my-gcs/bucket/chunks/b38/20cM/chunks_chrCNUMBER.txt",
            "rare_chunks": "gs://my-gcs/bucket/chunks/b38/4cM/chunks_chrCNUMBER.txt",
            "genetic_maps": "gs://my-gcs/bucket/maps/b38/chrCNUMBER.b38.gmap.gz"
        }


Next thing to do is update the :code:`nextflow.config` file. The only things you need to change are :code:`workDir` and
:code:`google.project`, and sometimes :code:`google.location`

    .. code-block:: sh

        workDir = 'gs://my-gcs/bucket/test_data/GWASpy/work'

        process {
          executor = 'google-batch'
          errorStrategy = { task.exitStatus==null ? 'retry' : 'terminate' }
          maxRetries = 3
        }

        profiles {
            gbatch {
              google.project = 'my-batch-billing-project'
              google.location = 'us-central1'
              batch.spot = true
            }
        }

Now you can easily run both phasing and imputation using the following command

    .. code-block:: sh

        ./nextflow run main.nf -c nextflow.config -profile gbatch -params-file params.json

5. Low-coverage WGS imputation using GLIMPSE
############################################

**COMING VERY SOON**
