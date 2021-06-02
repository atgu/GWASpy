# this is the script used to filter the reference 1KG+HGDP data used in PCA

import hail as hl

hl.init(default_reference='GRCh38')

ref_mt = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/hgdp_tgp_postQC.mt')

print("\nInitial number of SNPs before filtering: {}".format(ref_mt.count_rows()))
filtered_ref = hl.variant_qc(ref_mt)

print("filtering")

filtered_ref = filtered_ref.filter_rows((filtered_ref.variant_qc.AF[0] > 0.05) & (filtered_ref.variant_qc.AF[0] < 0.95))
print("\nNumber of SNPs after MAF filtering: {}".format(filtered_ref.count_rows()))

filtered_ref = filtered_ref.filter_rows(filtered_ref.variant_qc.call_rate > 0.999)
print("\nNumber of SNPs after Call Rate filtering: {}".format(filtered_ref.count_rows()))

# print("repartitioning")
# filtered_ref = filtered_ref.repartition(n_partitions=100, shuffle=True)

print("writing filtererd mt")
filtered_ref.write('gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_filtered_maf_5_GRCh38.mt', overwrite=True)
print("Done filtering")

print("getting sample information")

mt = hl.read_matrix_table('gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_filtered_maf_5_GRCh38.mt')

cols_ht = mt.cols()

pops = cols_ht.select(cols_ht.hgdp_tgp_meta.Study.region)

df = pops.to_pandas()

df.columns = ['Sample', 'SuperPop']

old_pops_labs = ['Africa', 'America', 'Central_South_Asia', 'East_Asia', 'Europe', 'Middle_East', 'Oceania', 'SAS']
new_pops_labs = ['AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID', 'OCE', 'CSA']
df['SuperPop'] = df['SuperPop'].replace(old_pops_labs, new_pops_labs)

print(df['SuperPop'].value_counts())

print("exporting sample metadata")
df.to_csv('gs://african-seq-data/hgdp_tgp/gwaspy_pca_ref/hgdp_1kg_sample_info.tsv', sep='\t', index=False)

