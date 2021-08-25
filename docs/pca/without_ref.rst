================================
Running PCA without a reference
================================

Put information on running PCA without a ref

Below is a code on how you can run PCA without a reference via the command-line or inside a Python script Use

#. Command line

    .. code-block:: sh

        pca --data-dirname data/ --data-basename 1kg_annotated --out-dir data/ --input-type hail --reference grch37

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="data/", data_basename="1kg_annotated",  out_dir="data/", input_type="hail", reference="GRCh37")