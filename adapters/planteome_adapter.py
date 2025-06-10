from functools import lru_cache
import hashlib
from enum import Enum, auto
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherNode, BioCypherEdge
from adapters.genome_adapter import GenomeAdapter
import re

logger.debug(f"Loading module {__name__}.")

class PlanteomeAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide
    """
    GENE_FUNCTION=auto()
    
    
class PlanteomeAdapterGeneFunctionField(Enum):
    """
    Define possible fields the adapter can provide for gene function
    """
    
    ID="Term"
    NAME="Name"
    DESCRIPTION="Description"
    EVIDENCE="Evidence"
    REFERENCE="Reference"
    CLASSEVIDENCE="ClassEvidence"
    
class PlanteomeAdapterEdgeType(Enum):
    """
    Enum for the types of the planteome adapter
    """
    
    FUNCTIONAL_ASSOCIATION="annotated_with"
    
class PlanteomeAdapter:
    """
    Planteome Adapter. Import, filter and create terms nodes and link them to genes

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
        self._set_types_and_fields(
            node_types,
            node_fields,
            edge_types,
        )
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the given file and extract terms.
        """
        logger.info("Preprocessing Planteome data.")

        # load data
        self.data = self._filter_input_planteome()
        
        # extract precursors (unique entities of `precursor` column)
        self.functions = self.data[["Term","Name","Description","Reference","Evidence","ClassEvidence"]].drop_duplicates()

        
    def get_nodes(self):
        """

        Returns a generator of BioCypher node objects for node types specified
        in the adapter constructor.

        Returns:
            Generator of BioCypher node objects.

        """

        logger.info("Generating nodes.")

        for _, row in self.functions.iterrows():
            node_id=row["Term"]
            name=row["Name"]
            description=row["Description"]
            reference=row["Reference"]
            evidence=row["Evidence"]
            classEvidence=row["ClassEvidence"]
            properties = {
                "preferred_id":node_id,
                "name":name,
                "description":description,
                "reference":reference,
                "evidence":evidence,
                "classEvidence":classEvidence
            }
        
            yield BioCypherNode(
                node_id=node_id,
                node_label="gene function",
                properties=properties,
            )
        
        
        
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.data.iterrows():
            # extract source and target
            source_id = row["OLN"]
            target_id = row["Term"]


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
                relationship_label="annotated_with",
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
            self.node_types = [type for type in PlanteomeAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in PlanteomeAdapterGeneFunctionField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in PlanteomeAdapterEdgeType]

    @lru_cache(maxsize=None)
    def _prefix_gene(self, string):
        return f"gene:{string}"
    
    def _extract_gene_id_planteome(self,row):
    # Check if gene ID is in Column2
        if pd.notna(row[2]):
            match = re.search(r'Solyc\d+g\d+', row[2])
            if match:
                return match.group()
        # Check if gene ID is in Column10
        if pd.notna(row[10]):
            match = re.search(r'Solyc\d+g\d+', row[10])
            if match:
                return match.group()
        #If 0 matches
        return None
    
    def _filter_input_planteome(self):
    
        Annotations=pd.read_csv('download/planteomeGene/select', header=None, sep='\t')
        Planteome_terms=pd.read_csv('download/planteomeTerms/select', header=None, sep='\t')
        
        #Set 'OLN' column index ont the column with the goog info
        Annotations['OLN'] = Annotations.apply(self._extract_gene_id_planteome, axis=1)
        
        Sly_genes=Annotations[[5,6,8,4,'OLN']]
        
        Associate_sly_genes_term = pd.merge(Sly_genes,Planteome_terms, left_on=4, right_on=0, how='inner')

        Associate_sly_genes_term.rename(columns={0:'Term',1:'Name',2:'Description',5:'Reference',6:'Evidence',8:'ClassEvidence'},inplace=True)

        Associate_sly_genes_term.drop(columns=4, inplace=True)

        Associate_sly_genes_term.drop(Associate_sly_genes_term[Associate_sly_genes_term['OLN']=='None'].index,inplace=True)
        Associate_sly_genes_term.dropna(inplace=True)
        
        Associate_sly_genes_term['Name'] = Associate_sly_genes_term['Name'].str.replace("'","''")
        Associate_sly_genes_term['Description'] = Associate_sly_genes_term['Description'].str.replace("'","''")
        Associate_sly_genes_term['OLN'] = Associate_sly_genes_term['OLN'].str.replace("'","''")
        Associate_sly_genes_term['Reference'] = Associate_sly_genes_term['Reference'].str.replace("'","''")
        Associate_sly_genes_term['Evidence'] = Associate_sly_genes_term['Evidence'].str.replace("'","''")
        Associate_sly_genes_term['ClassEvidence'] = Associate_sly_genes_term['ClassEvidence'].str.replace("'","''")
        
        Planteome_genes_sly_filt=GenomeAdapter(self.genome_path).filter_input_genome(Associate_sly_genes_term, 'OLN')

        
        return Planteome_genes_sly_filt