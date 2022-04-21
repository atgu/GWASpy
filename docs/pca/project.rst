================================
Project PCA (with a reference)
================================

You can leverage reference panel information to see how samples in your data cluster on a "global" scale.
PCs are computed using 1KG+HGDP dataset as a reference panel, and then samples in the input dataset are projected onto the 1KG+HGDP PC space.
A random forest classifier model, adopted from gnomAD, is then used to assign population ancestries in the input dataset

Below is a code on how you can run projection PCA via the command-line or inside a Python script Use

#. Command line

    .. code-block:: sh

        pca --data-dirname data/ --data-basename 1kg_annotated --out-dir data/ --input-type hail --reference grch37 --pca-type project

#. Python (inside a Python script)

    .. code-block:: python

        import gwaspy.pca as pca
        pca.pca.pca(data_dirname="data/", data_basename="1kg_annotated",  out_dir="data/",
                    input_type="hail", reference="GRCh37", pca_type="project")