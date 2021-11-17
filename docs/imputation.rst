.. _sec-imputation:

===================
Genotype Imputation
===================

Genotype imputation is a process of estimating missing genotypes from the haplotype or genotype reference panel. It
allows you to accurately evaluate the evidence for association at genetic markers that are not directly genotyped.
GWASpy has a module, :code:`imputation`, for running imputation using IMPUTE5. Because imputation is a computationally
intensive task, we run it on multiple chunks in parallel, then merge the imputed chunks together at the end. This is why
the module is divided into two parts: (1) :code:`impute`; (2) :code:`concat`. Below are examples of how to run imputation

1. Imputation
#############
#. Command line

    .. code-block:: sh

        imputation --input-vcfs gs://path/to/vcf_files.txt --samples-file gs://path/to/female_samples.txt --out-dir gs://path/to/output/dir --billing-project project-name --bucket bucket-associated-with-project --run impute --n-samples integer_number_of_samples

#. Python (inside a Python script)

    .. code-block:: python

            import gwaspy.imputation as impute
            impute.imputation.genotype_imputation(input_vcfs = 'gs://path/to/vcf_files.txt',
                      females_file: str = gs://path/to/female_samples.txt,
                      n_samples: int = integer_number_of_samples,
                      n_panel_samples: int = 4099,
                      buffer_region: int = 250,
                      local: bool = False,
                      billing_project = 'project-name',
                      bucket = 'bucket-associated-with-project',
                      memory: str = 'highmem',
                      cpu: int = 8,
                      run: str = 'impute'
                      output_type: str = 'bcf',
                      out_dir = 'gs://path/to/output/dir')

2. Concat
#########
After running imputation, we have to merge the imputed chunks together into full chromosomes. Here are examples below
#. Command line

    .. code-block:: sh

        imputation --input-vcfs gs://path/to/vcf_files.txt --samples-file gs://path/to/female_samples.txt --out-dir gs://path/to/output/dir --billing-project project-name --bucket bucket-associated-with-project --run concat --n-samples integer_number_of_samples

#. Python (inside a Python script)

    .. code-block:: python

            import gwaspy.imputation as impute
            impute.imputation.genotype_imputation(input_vcfs = 'gs://path/to/vcf_files.txt',
                      females_file: str = gs://path/to/female_samples.txt,
                      n_samples: int = integer_number_of_samples,
                      n_panel_samples: int = 4099,
                      buffer_region: int = 250,
                      local: bool = False,
                      billing_project = 'project-name',
                      bucket = 'bucket-associated-with-project',
                      memory: str = 'highmem',
                      cpu: int = 8,
                      run: str = 'concat'
                      output_type: str = 'bcf',
                      out_dir = 'gs://path/to/output/dir')



Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--input-vcfs`
     - Path to where text file containing VCF(s) for target genotypes paths is
   * - :code:`--samples-file`
     - Text file with list of FEMALE samples, one sample ID each line, that are in the dataset. This is crucial for chromosome X imputation as the data is split by sex
   * - :code:`--local`
     - Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud
   * - :code:`--billing-project`
     - Billing project to be used for the job(s)
   * - :code:`--bucket`
     - Bucket associated with the billing project
   * - :code:`--memory`
     - Memory to use for imputation. Options: [:code:`lowmem`, :code:`standard`, :code:`highmem`]. Default is :code:`highmem`
   * - :code:`--cpu-concat`
     - CPU to use for the concatenation step. Default is 8
   * - :code:`--n-samples`
     - Total number of samples in your dataset. We use this to estimate some of the job resources like storage.
   * - :code:`--buffer-region`
     - Buffer region to be used during imputation. This helps prevent imputation quality from deteriorating near the edges of the region. Default is 250 KB
   * - :code:`--run`
     - Process to run. Options: [:code:`impute`, :code:`concat`]. Default is :code:`impute`
   * - :code:`--out-type`
     - Output type. Options: [:code:`bcf`, :code:`vcf`]. Default is :code:`bcf` [HIGHLY RECOMMENDED SINCE BCFs ARE GENERALLY MORE EFFICIENT TO WORK WITH AND TAKE UP LESS SPACE]
   * - :code:`--out-dir`
     - Path to where output files will be saved

Output
######
The resulting output is a VCF file per chromosome with imputed genotypes.

.. note::
    Concatenating BCFs from imputation by chromosome is slower when the output is VCF compared to a BCF. The size may
    also differ significantly between BCF and VCF.