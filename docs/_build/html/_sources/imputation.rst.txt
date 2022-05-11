.. _sec-imputation:

===================
Genotype Imputation
===================

Genotype imputation is a process of estimating missing genotypes from the haplotype or genotype reference panel. It
allows you to accurately evaluate the evidence for association at genetic markers that are not directly genotyped.
GWASpy has a module, :code:`imputation`, for running imputation using IMPUTE5. Because imputation is a computationally
intensive task, we run it on multiple chunks in parallel, then merge the imputed chunks together at the end. This is why
the module is divided into two parts: (1) :code:`impute`; (2) :code:`concat`. Below are examples of how to run imputation

A. Run imputation+concat in a single command
#############################################
#. Command line

    .. code-block:: sh

        imputation --input-vcf gs://path/to/file.vcf.bgz --samples-file gs://path/to/female_samples.txt --out-dir gs://path/to/output/dir --billing-project project-name --run impute --n-samples integer_number_of_samples

#. Python (inside a Python script)

    .. code-block:: python

            import gwaspy.imputation as impute
            impute.imputation.genotype_imputation(input_vcfs = 'gs://path/to/file.vcf.bgz',
                      females_file: str = gs://path/to/female_samples.txt,
                      n_samples: int = integer_number_of_samples,
                      n_panel_samples: int = 4099,
                      buffer_region: int = 250,
                      local: bool = False,
                      billing_project = 'project-name',
                      memory: str = 'highmem',
                      cpu: int = 8,
                      stages: str = 'impute,concat'
                      output_type: str = 'bcf',
                      out_dir = 'gs://path/to/output/dir')

B. Run imputation and concat in separate commands
##################################################
If you want to run impute or concat as separate steps, you can set the :code:`--stages` (command-line)/:code:`stages` (Python script)
argument as impute or concat. It's important to note though that if you want to run things this way, the impute step should
always be run before concat as GWASpy uses results from the :code:`impute` stage for :code:`concat`



Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--input-vcf`
     - Path to where the VCF for target genotypes paths is
   * - :code:`--samples-file`
     - Text file with list of FEMALE samples, one sample ID each line, that are in the dataset. This is crucial for chromosome X imputation as the data is split by sex
   * - :code:`--local`
     - Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud
   * - :code:`--billing-project`
     - Billing project to be used for the job(s)
   * - :code:`--memory`
     - Memory to use for imputation. Options: [:code:`lowmem`, :code:`standard`, :code:`highmem`]. Default is :code:`highmem`
   * - :code:`--cpu-concat`
     - CPU to use for the concatenation step. Default is 8
   * - :code:`--n-samples`
     - Total number of samples in your dataset. We use this to estimate some of the job resources like storage.
   * - :code:`--buffer-region`
     - Buffer region to be used during imputation. This helps prevent imputation quality from deteriorating near the edges of the region. Default is 250 KB
   * - :code:`--stages`
     - Process to run. Default is :code:`impute,concat`
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