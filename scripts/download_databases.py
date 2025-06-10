from biocypher import BioCypher, FileDownload
import os
import tarfile
from biocypher._logger import logger
import subprocess


bc=BioCypher()

genome=FileDownload(
    name="genome",
    url_s=["https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/annotation/ITAG4.1_release/ITAG4.1_descriptions.txt",
           "https://solgenomics.net/ftp/genomes/Solanum_lycopersicum/annotation/ITAG4.1_release/ITAG4.1_gene_models.gff"],
    lifetime=1000
)


mirbase= FileDownload(
    name='mirbase',
    url_s=['https://www.mirbase.org/download/miRNA.dat',
           'https://www.mirbase.org/download/mature.fa',
           'https://www.mirbase.org/download/hairpin.fa'],
    lifetime=1000
)

plantregmap=FileDownload(
    name='plantregmap',
    url_s='https://plantregmap.gao-lab.org/download_ftp.php?filepath=08-download/Solanum_lycopersicum/binding/regulation_merged_Sly.txt',
    lifetime=1000
)

planttfdb=FileDownload(
    name='planttfdb',
    url_s='https://planttfdb.gao-lab.org/download/TF_list/Sly_TF_list.txt.gz',
    lifetime=1000
)

tardb= FileDownload(
    name='tardb',
    url_s='http://www.biosequencing.cn/TarDB/download/sly.zip',
    lifetime=1000
)

dpmind=FileDownload(
    name='dpmind',
    url_s='https://cbi.njau.edu.cn/DPMIND/sequence/verified_results/Solanum_lycopersicum.xlsx',
    lifetime=1000
)

string = FileDownload(
    name='string',
    url_s=['https://stringdb-downloads.org/download/protein.links.v12.0/4081.protein.links.v12.0.txt.gz',
           'https://stringdb-downloads.org/download/protein.aliases.v12.0/4081.protein.aliases.v12.0.txt.gz'],
    lifetime=1000
)


kegg = FileDownload(
    name='kegg',
    url_s=['https://rest.kegg.jp/link/sly/pathway',
           'https://rest.kegg.jp/list/pathway/sly',
           'https://rest.kegg.jp/conv/sly/uniprot'],
    lifetime=1000
)


#.. Planteome

userTaxon=('Solanum lycopersicum').split(' ')
taxon=userTaxon[0]+'%20'+userTaxon[1]

url='https://browser.planteome.org/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000000&start=0&fl=source,bioentity_internal_id,bioentity_label,qualifier,annotation_class,reference,evidence_type,evidence_with,aspect,bioentity_name,synonym,type,taxon,date,assigned_by,annotation_extension_class,bioentity_isoform&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=taxon_label:%22'
url=url+taxon+'%22&facet.field=source&facet.field=assigned_by&facet.field=aspect&facet.field=evidence_type_closure&facet.field=qualifier&facet.field=taxon_label&facet.field=type&facet.field=annotation_class_label&facet.field=regulates_closure_label&facet.field=annotation_extension_class_closure_label&q=*:*'


planteomeGene = FileDownload(
    name='planteomeGene',
    url_s=url,
    lifetime=1000
)

planteomeTerms = FileDownload(
    name='planteomeTerms',
    url_s='https://browser.planteome.org/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=1000000&start=0&fl=annotation_class,annotation_class_label,description&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class=%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22ontology_class%22&fq=is_obsolete:%22false%22&facet.field=source&facet.field=subset&facet.field=regulates_closure_label&facet.field=is_obsolete&q=*:*',
    lifetime=1000
)

def check_and_unzip(folder_name, tar_file_name):
    """Unzip the raw databases

    Args:
        folder_name (str): Folder name of the unzip file
        tar_file_name (str): Zip file
    """
    try:
        # Check if the folder exists and is not empty
        if os.path.exists(folder_name) and os.path.isdir(folder_name) and os.listdir(folder_name):
            print(f"The folder '{folder_name}' is already unzipped.")
            return True
        else:
            # Unzip the tar.gz file
            if tarfile.is_tarfile(tar_file_name):
                with tarfile.open(tar_file_name, "r:gz") as tar:
                    tar.extractall(folder_name)
                print(f"The folder '{folder_name}' has been unzipped.")
                return True
            else:
                print(f"The file '{tar_file_name}' is not a valid tar.gz file.")
                return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False

if __name__ == "__main__":
    folder_name = "database_downloads"
    tar_file_name = "database_downloads.tar.gz"
    check_and_unzip(folder_name, tar_file_name)
    bc.download(genome)
    bc.download(planteomeGene)
    bc.download(planteomeTerms)
    bc.download(plantregmap)
    bc.download(planttfdb)
    bc.download(tardb)
    bc.download(dpmind)
    bc.download(string)
    bc.download(mirbase)
    bc.download(kegg)
