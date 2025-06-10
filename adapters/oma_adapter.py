from functools import lru_cache
import hashlib
from enum import Enum, auto
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherNode, BioCypherEdge
from adapters.genome_adapter import GenomeAdapter
from omadb.OMARestAPI import Client

def _patch_pairwise_call():
        """
        Patches the `__call__` method of the `PairwiseRelations` class from the `omadb.OMARestAPI` module.
        This function replaces the original `__call__` method of the `PairwiseRelations` class with a new implementation
        that modifies the behavior of the method to make a customized API request.
        The new `__call__` method sends a request to retrieve pairwise relations including xrefs.
        
        The patched method simply add the following parameters to the original method:
            - `xrefs` (str, optional): Cross-references for filtering results. Defaults to None.
        Returns:
            None: This function modifies the `PairwiseRelations` class in place.
        Note:
            This patch should be used with caution as it modifies the behavior of the `PairwiseRelations` class globally.
        """
        from omadb.OMARestAPI import PairwiseRelations
        original_call = PairwiseRelations.__call__
        def __new_call__(self, genome_id1, genome_id2, chr1=None, chr2=None, xrefs=None,
                        rel_type=None, progress=False):
        
                return self._client._request(action=['pairs', genome_id2],
                                        subject=genome_id1,
                                        params={'chr1': chr1,
                                                'chr2': chr2,
                                                'rel_type': rel_type,
                                                'xrefs':xrefs},
                                        paginated=True,
                                        progress_desc=('Loading pairs'
                                                        if progress else None))
        
        PairwiseRelations.__call__ = __new_call__

_patch_pairwise_call()


logger.debug(f"Loading module {__name__}.")

class OmaAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide
    """
    ATHALIANA_GENE=auto()
    
    
class OmaAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for pathways
    """
    CANONICAL_ID="canonical_id"
    SPECIES="species"
    START="start"
    END="end"
    STRAND="strand"
    
class OmaAdapterEdgeType(Enum):
    """
    Enum for the types of the Oma adapter
    """
    
    GENE_TO_GENE_HOMOLOGY_ASSOCIATION="orthologous_to"
    
class OmaAdapterGeneToGeneHomologyAssociationEdgeField(Enum):
    """
    Define possible fields the adapter can provide for homology edges
    """
    
    RELATION_TYPE = "rel_type"
    DISTANCE = "distance"
    SCORE = "score"
    
class OmaAdapter:
    """
    Oma Adapter. Import, filter and create pathways nodes and link them to genes

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
        edge_fields: Optional[list] = None,
    ):
        self.genome_path = genome_path
        self._set_types_and_fields(
            node_types,
            node_fields,
            edge_types,
            edge_fields,
        )
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the API and extract the information.
        """
        
        self.data = self._filter_input_oma()
        self.ath_genes = self.data[['Ath_OLN', "Ath_canonicalid", 'species', 'start', 'end', 'strand']].drop_duplicates()
        self.homology = self.data[['OLN', 'Ath_OLN', 'rel_type', 'distance', 'score']].drop_duplicates()

        
    def get_nodes(self):
        """

        Returns a generator of BioCypher node objects for node types specified
        in the adapter constructor.

        Returns:
            Generator of BioCypher node objects.

        """

        logger.info("Generating nodes.")

        for _, row in self.ath_genes.iterrows():
            node_id=row["Ath_OLN"]
            canonical_id=row["Ath_canonicalid"]
            species=row["species"]
            start=row["start"]
            end=row["end"]
            strand=row["strand"]
            properties = {
                "name":node_id,
                "canonical_id":canonical_id,
                "species":species,
                "start":start,
                "end":end,
                "strand":strand,
            }
        
            yield BioCypherNode(
                node_id=node_id,
                node_label="Athaliana gene",
                properties=properties,
            )
        
        
        
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.homology.iterrows():
            # extract source and target
            source_id = row["OLN"]
            target_id = row["Ath_OLN"]
           
           #extract edge properties
            
            properties = {}
            
            if ( OmaAdapterGeneToGeneHomologyAssociationEdgeField.RELATION_TYPE in self.edge_fields ):
                properties["relation_type"] = row["rel_type"]
            
            if ( OmaAdapterGeneToGeneHomologyAssociationEdgeField.DISTANCE in self.edge_fields ):
                properties["distance"] = row["distance"]
                
            if ( OmaAdapterGeneToGeneHomologyAssociationEdgeField.SCORE in self.edge_fields ):
                properties["score"] = row["score"]
            

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
                relationship_label="orthologous_to",
                properties=properties,
            )
       
    def _set_types_and_fields(
        self,
        node_types,
        node_fields,
        edge_types,
        edge_fields,
    ):
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in OmaAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in OmaAdapterGeneField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in OmaAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in OmaAdapterGeneToGeneHomologyAssociationEdgeField]

    @lru_cache(maxsize=None)
    def _prefix_gene(self, string):
        return f"gene:{string}"
    
    def _filter_input_oma(self):
        """ Filter input oma gene to match the Genome

        Returns:
            oma_filtered (pd.DataFrame): Dataframe with pathway genes association after genome filtering of the genes
        """
    
        oma = self.combine_oma(self.gene_in_pathway_file, self.gene_to_uniprot_file)
        
        # oma_genes['OLN']=oma_genes['OLN'].str.split('.').str.get(0)
        
        # print('oma interactions:', oma[['OLN','Pathway_ID']].drop_duplicates().shape[0])
        
        oma_filtered=GenomeAdapter(self.genome_path).filter_input_genome(oma, 'OLN')
        
        # print('oma interactions after filtering:', oma_filtered[['OLN','Pathway_ID']].drop_duplicates().shape[0])
        
        return oma_filtered
    
    def _get_oma_information(self):
        
        logger.info("Retrieving OMA information.")
        
        c = Client()
        # get the pairwise relations between Sly and Ath
        tomato_arath_pairs = c.pairwise('SOLLC','ARATH', progress=True, xrefs='SourceAC').as_dataframe()
        
        return tomato_arath_pairs
    
    def _parse_oma(self):
        
        oma_sly_ath = self._get_oma_information()
        #Parse the information
        # Extract and flatten first
        entry1_df = pd.json_normalize(oma_sly_ath['entry_1'])
        entry2_df = pd.json_normalize(oma_sly_ath['entry_2'])

        # Optionally add prefixes
        entry1_df = entry1_df.add_prefix('Sly_')
        entry2_df = entry2_df.add_prefix('Ath_')

        # Drop the nested columns
        oma_sly_ath_parsed = oma_sly_ath.drop(columns=['entry_1', 'entry_2'])

        # Reset all indexes to ensure alignment
        oma_sly_ath_parsed = oma_sly_ath_parsed.reset_index(drop=True)
        entry1_df = entry1_df.reset_index(drop=True)
        entry2_df = entry2_df.reset_index(drop=True)

        # Concatenate safely
        oma_sly_ath_parsed = pd.concat([entry1_df, entry2_df, oma_sly_ath_parsed], axis=1)
    
        
        return oma_sly_ath_parsed
    
    def _retrieve_info_from_parsed(self):
        
        oma_parsed = self._parse_oma()
        # Extract the relevant columns
        columns = ["Sly_xrefs.SourceAC", "Ath_canonicalid","Ath_xrefs.SourceAC", "Ath_locus.start", "Ath_locus.end",
                   "Ath_locus.strand","Ath_species.species", "rel_type", "distance", "score"]
        
        # Select the relevant columns
        oma_parsed = oma_parsed[columns]
        # Rename the columns
        oma_parsed = oma_parsed.rename(columns={
            "Sly_xrefs.SourceAC": "OLN",
            "Ath_xrefs.SourceAC": "Ath_OLN",
            "Ath_species.species": "species",
            "Ath_locus.start": "start",
            "Ath_locus.end": "end",
            "Ath_locus.strand": "strand",
            })
        
        oma_parsed['OLN']=oma_parsed['OLN'].str.split('.').str.get(0)
        
        return oma_parsed
    
    def _filter_input_oma(self):
        """ Filter input oma gene to match the Genome

        Returns:
            oma_filtered (pd.DataFrame): Dataframe with pathway genes association after genome filtering of the genes
        """
    
        oma = self._retrieve_info_from_parsed()
        
        # print('oma interactions:', oma[['OLN','Pathway_ID']].drop_duplicates().shape[0])
        
        oma_filtered=GenomeAdapter(self.genome_path).filter_input_genome(oma, 'OLN')
        
        # print('oma interactions after filtering:', oma_filtered[['OLN','Pathway_ID']].drop_duplicates().shape[0])
        
        return oma_filtered 
    
    
        
    
    