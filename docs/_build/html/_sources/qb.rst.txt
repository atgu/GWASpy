.. _sec-qb:

====================
Hail Query and Batch
====================

The four GWASpy modules use two different backends: :code:`preimp_qc` and :code:`pca` use Hail Query, while
:code:`phasing` and :code:`imputation` modules use Batch (Hail Batch for Broad users and nextflow for non-Broad users).
Hail Query is well-suited for manipulating large genomics data in a highly parallelised environments such as Dataproc.
`Batch <https://cloud.google.com/batch/docs/get-started>`_, on the other hand, is good for batch processing (scheduling,
queueing, and executing) workloads on Google Cloud resources.

All the instructions below assume the user has a Google account and an active (Google) Cloud billing account

Query
#####

For running the :code:`preimp_qc` and :code:`pca` modules, you need to start a Dataproc cluster. Hail has a command-line
tool, `hailctl <https://hail.is/docs/0.2/cloud/google_cloud.html>`_, for doing this and it is installed automatically when
you install Hail. We highly recommend setting a maximum age for the cluster (:code:`--max-age`), this will ensure the cluster is
automatically deleted after the specified time.

Below is how you can start a cluster with GWASpy pre-installed:

    .. code-block:: sh

       hailctl dataproc start my-cluster-name -region=us-central1 --packages gwaspy --max-age 4h

To shut down the cluster, you can run:

    .. code-block:: sh

        hailctl dataproc stop my-cluster-name --region=us-central1

Batch
#####

The :code:`phasing` and :code:`imputation` modules use Batch as the backend. For Broad users with a Hail Batch account,
there is no setup needed, you can proceed to running the modules. For non-Broad users, we have a nextflow implementation
of the modules that requires nextflow setup first. Follow the steps here to: `(1) install nextflow <https://www.nextflow.io/docs/latest/install.html#install-page>`_; and
`(2) setup Google Cloud Batch for nextflow <https://www.nextflow.io/docs/latest/google.html>`_
