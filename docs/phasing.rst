.. _sec-phasing:

=================
Haplotype Phasing
=================

Knowing the phase of a haplotype can allow us to impute low frequency variants, this makes haplotype phasing an
important step before genotype imputation. GWASpy has a module, :code:`phasing`, for performing phasing. Phasing can
be run with or without a reference panel using SHAPEIT5

GWASpy can handle both array and WGS data. For array data, the user can pass a VCF/BCF file with all the chromosomes,
then GWASpy will use SHAPEIT5 to phase the chromosomes in parallel. Since WGS has more variants, phasing will be parallelized across
multiple chunks in each chromosome. It's also important to note that phasing of WGS data includes phasing common
variants first, followed by phasing rare variants.

Another important aspect of phasing is the use of a reference panel. In many cases (small sample size), including a reference panel when
phasing improves accuracy. By default, GWASpy runs phasing without a reference panel, but there is an option to use a
reference panel as shown below.

Examples
########

**1. Without a reference panel**

    .. code-block:: sh

        phasing --input-vcf gs://path/to/file.vcf.bgz --output-filename outfilename.phased --out-dir gs://path/to/output/dir --genome-build GRCh38 --billing-project my-billing-project

**2. HGDP+1KG reference panel**

Set :code:`--vcf-ref` to  :code:`hgdp1kgp`

    .. code-block:: sh

        phasing --input-vcf gs://path/to/file.vcf.bgz --output-filename my_outfilename --out-dir gs://path/to/output/dir --genome-build GRCh38 --billing-project my-billing-project --vcf-ref hgdp1kgp

**3. Own reference panel**

.. note::
    1. If you're using your own reference panel, make sure the files are bgzip compressed.
    2. Chromosome X reference file must be named X and not 23

Say you have your reference panel files for each chromosomes stored in gs://ref_panel/ALL.chr{1..22,X}.vcf,
you would pass the path to :code:`--vcf-ref` as gs://ref_panel/ALL.chr\ **CNUMBER**\ .vcf.
GWASpy uses **CNUMBER** as a placeholder for the chromosomes. Then you can run phasing as:

    .. code-block:: sh

        phasing --input-vcf gs://path/to/file.vcf.bgz --output-filename outfilename.phased --out-dir gs://path/to/output/dir --genome-build GRCh38 --billing-project my-billing-project --vcf-ref gs://ref_panel/ALL.chrCNUMBER.vcf

.. note::
    For nextflow users, the idea is the same. The only difference is you have to update the params.json file. Examples
    are provided in the tutorial section of the documentation

Arguments and options
#####################

.. list-table::
   :widths: 15 50
   :header-rows: 1

   * - Argument
     - Description
   * - :code:`--input-vcf`
     - Path to where VCF file to be phased is
   * - :code:`--vcf-ref`
     - VCF file for reference haplotypes if phasing with a reference panel
   * - :code:`--pedigree`
     - Pedigree (PLINK FAM) file
   * - :code:`--local`
     - Type of service. Default is Service backend where jobs are executed on a multi-tenant compute cluster in Google Cloud
   * - :code:`--billing-project`
     - Billing project to be used for the job(s)
   * - :code:`--genome-build`
     - Genome reference build. Default is GRCh38. Options: [:code:`GRCh37`, :code:`GRCh38`]
   * - :code:`--data-type`
     - Array or WGS data. Default is array. Options: [:code:`array`, :code:`wgs`].
   * - :code:`--fill-tags`
     - Whether or not to add AC tag required by SHAPEIT5. Including :code:`--fill-tags`, in your command will enable this step
   * - :code:`--software`
     - Software to use for phasing. Options: [:code:`beagle`, :code:`shapeit`]. Default is :code:`shapeit`
   * - :code:`--output-filename`
     - Output filename without file extension
   * - :code:`--out-dir`
     - Path to where output files will be saved

Output
######
The resulting output is a VCF file per chromosome with phased haplotypes.
