================================
Normal PCA (without a reference)
================================

GWASpy allows you to run normal PCA without any reference panel

Below is a code on how you can run normal PCA without a reference via the command-line or inside a Python script Use

#. Command line

    .. code-block:: sh

        pca --data-dirname data/ --data-basename 1kg_annotated --out-dir data/ --input-type hail --reference grch37 --pca-type normal

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="data/", data_basename="1kg_annotated",  out_dir="data/",
                    input_type="hail", reference="GRCh37", pca_type="normal")