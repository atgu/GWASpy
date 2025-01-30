================================
Normal PCA (without a reference)
================================

GWASpy allows you to run normal PCA without any reference panel

Below is a code on how you can run normal PCA without a reference via the command-line or inside a Python script Use

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="gs://my-gcs/bucket/test_data/", data_basename="my_data_basename",
                    out_dir="gs://my-gcs/bucket/test_data/", input_type="my_input_type", reference="GRCh37",
                    pca_type="normal")

#. Command line

    .. code-block:: sh

        pca --data-dirname gs://my-gcs/bucket/test_data/ --data-basename my_data_basename --out-dir gs://my-gcs/bucket/test_data/--input-type my_input_type --reference grch37 --pca-type normal