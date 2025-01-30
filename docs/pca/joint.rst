================================
Joint PCA (with a reference)
================================

The joint PCA method works by first merging (joining), by locus and allele(s), the input dataset with the reference dataset.
This is followed by "normal" PCA on the merged dataset

Below is a code on how you can run joint PCA via the command-line or inside a Python script Use

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="gs://my-gcs/bucket/test_data/", data_basename="my_data_basename",
                    out_dir="gs://my-gcs/bucket/test_data/", input_type="my_input_type", reference="GRCh37",
                    pca_type="joint")

#. Command line

    .. code-block:: sh

        pca --data-dirname gs://my-gcs/bucket/test_data/ --data-basename my_data_basename --out-dir gs://my-gcs/bucket/test_data/--input-type my_input_type --reference grch37 --pca-type joint