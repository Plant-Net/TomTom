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

class MercatorAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide
    """
    PATHWAY=auto()
    
    
class MercatorAdapterGeneFunctionField(Enum):
    """
    Define possible fields the adapter can provide for pathways
    """
    
    ID="Pathway"
    
class MercatorAdapterEdgeType(Enum):
    """
    Enum for the types of the Mercator adapter
    """
    
    GENE_TO_PATHWAY_ASSOCIATION="involved_in"
    
class MercatorAdapter:
    """
    Mercator Adapter. Import, filter and create pathways nodes and link them to genes

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
        self.pathway_file = 'database_downloads/database_downloads/mercator/Mercator_annotation_Sly_4_1.txt'
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
    
        self.annotation = self._filter_input_Mercator()
        
        # extract precursors (unique entities of `precursor` column)
        self.pathways = self.annotation[["Pathway","Annotation"]].drop_duplicates()
        self.genes = self.annotation[["Pathway","OLN"]].drop_duplicates()

        
    def get_nodes(self):
        """

        Returns a generator of BioCypher node objects for node types specified
        in the adapter constructor.

        Returns:
            Generator of BioCypher node objects.

        """

        logger.info("Generating nodes.")

        for _, row in self.pathways.iterrows():
            node_id=row["Pathway"]
            properties = {
                "name":node_id
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
            target_id = row["Pathway"]


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
            self.node_types = [type for type in MercatorAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in MercatorAdapterGeneFunctionField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in MercatorAdapterEdgeType]

    @lru_cache(maxsize=None)
    def _prefix_gene(self, string):
        return f"gene:{string}"
    
    def _filter_input_Mercator(self):
        """ Filter input Mercator gene to match the Genome

        Returns:
            Mercator_filtered (pd.DataFrame): Dataframe with pathway genes association after genome filtering of the genes
        """
    
        mercator = self.read_pathways()
        
        
        mercator_filtered=GenomeAdapter(self.genome_path).filter_input_genome(mercator, 'OLN')
        
        
        return mercator_filtered
    

    def read_pathways(self):
        """Read and parse the pathways in the S. lycopersicum species from the file


        Returns:
            mercator (pd.DataFrame): DataFrame of the pathways Name and Annotation description.
        """
        
        mercator = pd.read_csv(self.pathway_file, sep="\t", quotechar="'")
        mercator["IDENTIFIER"]= mercator["IDENTIFIER"].str.split(".").str.get(0)
        mercator["IDENTIFIER"]= mercator["IDENTIFIER"].str.capitalize()
        mercator = mercator[["NAME","IDENTIFIER"]]
        mercator["Pathway"] = mercator["NAME"].str.split(".").str.get(0)
        mercator.rename(columns={"NAME":"Annotation", "IDENTIFIER": "OLN"}, inplace=True)
        mercator = mercator[mercator["Pathway"] != 'not assigned']
        
        return mercator
        
 

        
 
        
    
    