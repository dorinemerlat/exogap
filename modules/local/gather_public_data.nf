process GATHER_PUBLIC_DATA {

    input:
    path(proteins_count)
    path(proteins_from_proteomes_count)
    path(transcriptomes_count)

    output:
    path("infos_about_taxons.csv")

    script:
    """
    cat *proteins_count.tsv |sort > all_proteins.csv
    cat *proteomes_count.tsv |sort > all_proteins_from_proteoms.csv
    cat *transcriptomes_count.tsv |sort > all_transcriptomes.csv

    cut -f 3 -d ',' all_proteins_from_proteoms.csv > all_proteins_from_proteoms.csv.temp
    cut -f 3 -d ',' all_transcriptomes.csv > all_transcriptomes.csv.temp

    paste -d ',' all_proteins.csv all_proteins_from_proteoms.csv.temp all_transcriptomes.csv.temp > infos_about_taxons.csv

    # add header to infos_about_taxons.csv
    sed -i '1 i\\taxon_taxid,taxon_name,proteins_count,proteins_from_proteomes_count,transcriptomes_count' infos_about_taxons.csv
    """
}
