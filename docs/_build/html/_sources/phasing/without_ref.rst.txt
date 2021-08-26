====================================
Running Phasing without a reference
====================================

Users can run phasing without a reference panel using the code:`phasing` module

Below is a code on how you can run phasing without a reference panel either: (1) via the command-line;
or (2) inside a Python script

#. Command line

    .. code-block:: sh

        phasing --input-vcfs gs://path/to/vcf_files.txt --out-dir gs://path/to/output/directory --billing-project project-name --bucket bucket-associated-with-project

#. Python (inside a Python script)

    .. code-block:: python

            import gwaspy.phasing as phase
            phase.phasing.haplotype_phasing(input_vcfs = 'gs://path/to/vcf_files.txt',
                      software = 'shapeit',
                      reference= 'GRCh37',
                      out_dir = 'gs://path/to/output/dir',
                      billing_project = 'project-name',
                      bucket = 'bucket-associated-with-project')


