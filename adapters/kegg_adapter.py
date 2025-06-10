from functools import lru_cache
import hashlib
from enum import Enum, auto
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherNode, BioCypherEdge
from adapters.genome_adapter import GenomeAdapter
import os
import subprocess
from scripts.Uniprot_ID_mapping import *

logger.debug(f"Loading module {__name__}.")

class KeggAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide
    """
    PATHWAY=auto()
    
    
class KeggAdapterGeneFunctionField(Enum):
    """
    Define possible fields the adapter can provide for pathways
    """
    
    ID="Pathway_ID"
    DESCRIPTION="Pathway_Description"
    SPECIES="Species"
    
class KeggAdapterEdgeType(Enum):
    """
    Enum for the types of the Kegg adapter
    """
    
    GENE_TO_PATHWAY_ASSOCIATION="involved_in"
    
class KeggAdapter:
    """
    KEGG Adapter. Import, filter and create pathways nodes and link them to genes

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
    """

    def __init__(
        self,
        genome_path: str,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
    ):
        self.genome_path = genome_path
        self.pathway_file = 'download/kegg/sly'
        self.gene_in_pathway_file = 'download/kegg/pathway'
        self.gene_to_uniprot_file = 'download/kegg/uniprot'
        self._set_types_and_fields(
            node_types,
            node_fields,
            edge_types,
        )
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the given file and extract pathways and gene in pathways.
        """
        
        self.pathways = self.read_pathways(self.pathway_file)
        self.genes = self._filter_input_kegg()
        
        # extract precursors (unique entities of `precursor` column)
        self.pathways = self.pathways[["Pathway_Description","Pathway_ID","Species"]].drop_duplicates()
        self.genes = self.genes[["Pathway_ID","OLN"]].drop_duplicates()

        
    def get_nodes(self):
        """

        Returns a generator of BioCypher node objects for node types specified
        in the adapter constructor.

        Returns:
            Generator of BioCypher node objects.

        """

        logger.info("Generating nodes.")

        for _, row in self.pathways.iterrows():
            node_id=row["Pathway_ID"]
            name=row["Pathway_Description"]
            species=row["Species"]
            properties = {
                "preferred_id":node_id,
                "name":name,
                "species":species
            }
        
            yield BioCypherNode(
                node_id=node_id,
                node_label="pathway",
                properties=properties,
            )
        
        
        
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.genes.iterrows():
            # extract source and target
            source_id = row["OLN"]
            target_id = row["Pathway_ID"]


            # generate relationship id
            md5 = hashlib.md5(
                "".join(
                    [str(source_id), str(target_id)]
                ).encode("utf-8")
            ).hexdigest()

            # generate edge
            yield BioCypherEdge(
                relationship_id=md5,
                source_id=self._prefix_gene(source_id),
                target_id=target_id,
                relationship_label="involved_in",
            )
       
    def _set_types_and_fields(
        self,
        node_types,
        node_fields,
        edge_types,
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in KeggAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in KeggAdapterGeneFunctionField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in KeggAdapterEdgeType]

    @lru_cache(maxsize=None)
    def _prefix_gene(self, string):
        return f"gene:{string}"
    
    def _filter_input_kegg(self):
        """ Filter input kegg gene to match the Genome

        Returns:
            kegg_filtered (pd.DataFrame): Dataframe with pathway genes association after genome filtering of the genes
        """
    
        kegg = self.combine_kegg(self.gene_in_pathway_file, self.gene_to_uniprot_file)
        
        # kegg_genes['OLN']=kegg_genes['OLN'].str.split('.').str.get(0)
        
        # print('KEGG interactions:', kegg[['OLN','Pathway_ID']].drop_duplicates().shape[0])
        
        kegg_filtered=GenomeAdapter(self.genome_path).filter_input_genome(kegg, 'OLN')
        
        # print('KEGG interactions after filtering:', kegg_filtered[['OLN','Pathway_ID']].drop_duplicates().shape[0])
        
        return kegg_filtered
    
    def combine_kegg(self,gene_to_pathways_file,gene_to_uniprot_file):
        """Combine the KEGG information and convert the uniprot id to OLN id.

        Args:
            gene_to_pathways_file (str): File path of the KEGG genes associated to pathways obtained from the download script
            gene_to_uniprot_file (str): File path of the KEGG conversion of gene id to uniprot id obtained from the download script

        Returns:
            KEGG (pd.DataFrame): DataFrame of 4 column representing the KEGG information: gene_id Pathway_ID uniprot_id and OLN
        """

        gene_to_pathways = self.read_gene_in_pathway(gene_to_pathways_file)
        gene_to_uniprot = self.read_gene_uniprot(gene_to_uniprot_file)
        
        KEGG = pd.merge(gene_to_pathways, gene_to_uniprot, on='gene_id')
        
        list_uniprot = KEGG['uniprot_id'].unique().tolist()
        
        uniprot_mapping = self.map_uniprot_to_oln(list_uniprot)
        
        KEGG = pd.merge(KEGG,uniprot_mapping, left_on='uniprot_id', right_on='From')
        KEGG.drop(['From','To'],axis=1,inplace=True)
        
        return KEGG

    def read_pathways(self,list_pathways_file):
        """Read and parse the pathways in the S. lycopersicum species from the file

        Args:
            list_pathways_file (str): File path of the KEGG pathways for the Sly species obtained from the download script

        Returns:
            list_pathways (pd.DataFrame): DataFrame of the pathways ID, pathways description (name) and species
        """
        list_pathways = pd.read_csv(list_pathways_file, sep='\t', header=None)
        list_pathways.rename(columns={0:'Pathway_ID',1:'description'}, inplace=True)
        list_pathways[['Pathway_Description', 'Species']] = list_pathways['description'].str.extract(r'(.+?)\s-\s([A-Z][a-z]+.*)')
        list_pathways.drop('description',axis=1,inplace=True)
        
        return list_pathways
        
    def read_gene_in_pathway(self,gene_to_pathways_file):
        """Read and parse the gene in the KEGG pathways

        Args:
            gene_to_pathways_file (str): File path of the KEGG genes associated to pathways obtained from the download script

        Returns:
            gene_to_pathways (pd.DataFrame): DataFrame of the gene (ncbi_id) pathways association
        """
        gene_to_pathways = pd.read_csv(gene_to_pathways_file, sep='\t', header=None)
        gene_to_pathways.rename(columns={0:'Pathway_ID',1:'gene_id'}, inplace=True)
        gene_to_pathways['gene_id'] = gene_to_pathways['gene_id'].str.split(':').str.get(1)
        gene_to_pathways['Pathway_ID'] = gene_to_pathways['Pathway_ID'].str.split(':').str.get(1)

        return gene_to_pathways

    def read_gene_uniprot(self,gene_to_uniprot_file):
        """Read and parse the conversion from gene_id (ncbi_id) to uniprot_id

        Args:
            gene_to_uniprot_file (str): File path of the KEGG conversion of gene id to uniprot id obtained from the download script

        Returns:
            gene_to_uniprot (pd.Dataframe): DataFrame of the conversion of gene_id to uniprot_id in KEGG
        """
        gene_to_uniprot = pd.read_csv(gene_to_uniprot_file, sep='\t', header=None)
        gene_to_uniprot.rename(columns={0:'uniprot_id',1:'gene_id'},inplace=True)
        gene_to_uniprot['gene_id'] = gene_to_uniprot['gene_id'].str.split(':').str.get(1)
        gene_to_uniprot['uniprot_id'] = gene_to_uniprot['uniprot_id'].str.split(':').str.get(1)
        
        return gene_to_uniprot
        
    def map_uniprot_to_oln(self,uniprot_ids):
        """Perform the ID mapping from UniProt to OLN using the Uniprot API.
        
        Args:
            uniprot_ids (list): List of UniProt IDs to convert to OLN IDs.
        
        Returns:
            uniprot_mapping (pd.DataFrame): DataFrame of the conversion from UniProt IDs to OLN IDs.
        """
        # Example of how to use UniProt's API to convert UniProt IDs to OLN IDs.
        job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="Ensembl_Genomes", ids=uniprot_ids)
        
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results_map = get_id_mapping_results_search(link + "?format=tsv")
        
        # Process results into a mapping
        uniprot_mapping = get_data_frame_from_tsv_results(results_map)
        uniprot_mapping['OLN'] = uniprot_mapping['To'].str.split('.').str.get(0)
        
        return uniprot_mapping
        
    
    