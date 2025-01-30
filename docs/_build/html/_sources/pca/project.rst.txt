================================
Project PCA (with a reference)
================================

You can leverage reference panel information to see how samples in your data cluster on a "global" scale.
PCs are computed using 1KG+HGDP dataset as a reference panel, and then samples in the input dataset are projected onto the 1KG+HGDP PC space.
A random forest classifier model, adopted from gnomAD, is then used to assign population ancestries in the input dataset

Below is a code on how you can run projection PCA via the command-line or inside a Python script Use

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="gs://my-gcs/bucket/test_data/", data_basename="my_data_basename",
                    out_dir="gs://my-gcs/bucket/test_data/", input_type="my_input_type", reference="GRCh37",
                    pca_type="project")

#. Command line

    .. code-block:: sh

        pca --data-dirname gs://my-gcs/bucket/test_data/ --data-basename my_data_basename --out-dir gs://my-gcs/bucket/test_data/--input-type my_input_type --reference grch37 --pca-type project