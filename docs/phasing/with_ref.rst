=================================
Running Phasing with a reference
=================================

Put information on when including a reference panel is beneficial

Below is a code on how you can run phasing with a reference panel.

.. code-block:: sh

   phasing --input-vcfs gs://path/to/vcf_files.txt --vcf-ref gs://path/to/reference_panel.vcf --out-dir gs://path/to/output/directory
