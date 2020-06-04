##########################################
#  Ovine Lymph Node Liver Fluke RNA-seq  #
##########################################
 
# Author: Amalia Naranjo
# Accessed: 04/06/2020
# https://sparql.orthodb.org v10.1


# Get all sheep gene IDs
prefix : <http://purl.orthodb.org/>
select *
where {
?gene_sheep a :Gene; :name ?gene_sheep_name; :description ?sheep_description; up:organism/a [up:scientificName "Ovis aries"]; :xref [a :Xref; :xrefResource ?xref_sheep].
}

# Get all human gene IDs
prefix : <http://purl.orthodb.org/>
select *
where {
?gene_human a :Gene; :name ?gene_human_name; :description ?human_description; up:organism/a [up:scientificName "Homo sapiens"]; :xref [a :Xref; :xrefResource ?xref_human].
}

# Get sheep to human orthologs
prefix : <http://purl.orthodb.org/>
select ?og ?og_description ?gene_sheep ?gene_sheep_name ?gene_human ?gene_human_name
where {
?gene_sheep a :Gene.
?gene_human a :Gene.
?gene_sheep up:organism/a [up:scientificName "Ovis aries"].
?gene_human up:organism/a taxon:9606.
?gene_sheep :name ?gene_sheep_name.
?gene_human :name ?gene_human_name.
?gene_sheep :memberOf ?og.
?gene_human :memberOf ?og.
?og a :OrthoGroup; :name ?og_description.
}