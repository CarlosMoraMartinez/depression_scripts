
# To include in NF pipeline

humann_regroup_table -i humann3_merged_genetables.tsv -o humann3_merged_genetables_GO.tsv -g uniref90_go
humann_regroup_table -i humann3_merged_genetables.tsv -o humann3_merged_genetables_KO.tsv -g uniref90_ko
humann_regroup_table -i humann3_merged_genetables.tsv -o humann3_merged_genetables_RXN.tsv -g uniref90_rxn
humann_regroup_table -i humann3_merged_genetables.tsv -o humann3_merged_genetables_CAZY.tsv -c /home/ysanz/ddbb/humann_dbs/utility_mapping/utility_mapping/map_cazy_uniref90.txt.gz

humann_renorm_table --input humann3_merged_genetables_GO.tsv --output humann3_merged_genetables_GO_CPM.tsv --units cpm --update-snames
humann_renorm_table --input humann3_merged_genetables_KO.tsv --output humann3_merged_genetables_KO_CPM.tsv --units cpm --update-snames
humann_renorm_table --input humann3_merged_genetables_RXN.tsv --output humann3_merged_genetables_RXN_CPM.tsv --units cpm --update-snames
humann_renorm_table --input humann3_merged_genetables_CAZY.tsv --output humann3_merged_genetables_CAZY_CPM.tsv --units cpm --update-snames
humann_renorm_table --input humann3_merged_abundances.tsv --output humann3_merged_abundances_CPM.tsv --units cpm --update-snames

humann_rename_table --input humann3_merged_genetables_RXN_CPM.tsv --output humann3_merged_genetables_RXN_CPM_namedPWY.tsv --names metacyc-pwy
humann_rename_table --input humann3_merged_genetables_KO_CPM.tsv --output humann3_merged_genetables_KO_CPM_named.tsv --names kegg-pathway
humann_rename_table --input humann3_merged_genetables_GO_CPM.tsv --output humann3_merged_genetables_GO_CPM_named.tsv --names go
humann_rename_table --input humann3_merged_abundances_CPM.tsv --output humann3_merged_abundances_CPM_namedPWY.tsv --names metacyc-pwy

humann_rename_table --input humann3_merged_genetables_RXN_CPM.tsv --output humann3_merged_genetables_RXN_CPM_namedRXN.tsv --names metacyc-rxn
humann_rename_table --input humann3_merged_abundances_CPM.tsv --output humann3_merged_abundances_CPM_namedRXN.tsv --names metacyc-rxn





