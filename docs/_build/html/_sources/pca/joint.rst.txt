================================
Joint PCA (with a reference)
================================

The joint PCA method works by first merging (joining), by locus and alleles, the input dataset with the reference dataset.
This is followed by "normal" PCA on the merged dataset

Below is a code on how you can run joint PCA via the command-line or inside a Python script Use

#. Command line

    .. code-block:: sh

        pca --data-dirname data/ --data-basename 1kg_annotated --out-dir data/ --input-type hail --reference grch37 --pca-type joint

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="data/", data_basename="1kg_annotated",  out_dir="data/",
                    input_type="hail", reference="GRCh37", pca_type="joint")