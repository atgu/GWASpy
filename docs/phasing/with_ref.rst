=================================
Running Phasing with a reference
=================================

Using a reference panel for phasing can improve the accuracy of phasing. The code:`phasing` module can be run with a
reference panel

Below is a code on how you can run phasing with a reference panel via:

#. Command line

    .. code-block:: sh

        phasing --input-vcfs gs://path/to/vcf_files.txt --vcf-ref gs://path/to/reference_panel.vcf --out-dir gs://path/to/output/dir --billing-project project-name --bucket bucket-associated-with-project

#. Python (inside a Python script)

    .. code-block:: python

            import gwaspy.phasing as phase
            phase.phasing.haplotype_phasing(input_vcfs = 'gs://path/to/vcf_files.txt',
                      vcf_ref = 'gs://path/to/reference_panel.vcf'
                      software = 'shapeit',
                      reference= 'GRCh37',
                      out_dir = 'gs://path/to/output/dir',
                      billing_project = 'project-name',
                      bucket = 'bucket-associated-with-project')
