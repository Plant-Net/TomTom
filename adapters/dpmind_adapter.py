from functools import lru_cache
import hashlib
from enum import Enum
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherEdge
from adapters.genome_adapter import GenomeAdapter
from adapters.miRBase_adapter import MirbaseAdapter

logger.debug(f"Loading module {__name__}.")
    
class DpmindAdapterEdgeType(Enum):
    """
    Enum for the types of the MTI adapter
    """
    
    MICRORNA_INHIBITION="represses_dpmind"
    
class DpmindAdapterMicrornaInhibitionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for microrna inhibition edges
    """
    
    DEGRADOME_SOURCE_DPMIND = "degradome"
    CATEGORY_DPMIND = "category"
    
class DpmindAdapter:
    """
    DPMIND adapter. Import, filter and create links between microRNAs and targets.

    Args:
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
        file_path: str,
        genome_path: str,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self.file_path = file_path
        self.genome_path = genome_path
        self._set_types_and_fields(
            edge_types,
            edge_fields,
        )
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the given file and extract micrornas and target mrnas.
        """
        logger.info("Preprocessing DPMIND data.")

        # load data
        self.data = self._filter_input_dpmind()
        
            
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.data.iterrows():
            # extract source and target
            source_id = row["mir_name"]
            target_id = row["tar_OLN"]
            
            #extract edge properties
            
            properties = {}
            
            if ( DpmindAdapterMicrornaInhibitionEdgeField.DEGRADOME_SOURCE_DPMIND in self.edge_fields ):
                properties["degradome source dpmind"] = row["degradome"]
            
            if ( DpmindAdapterMicrornaInhibitionEdgeField.CATEGORY_DPMIND in self.edge_fields ):
                properties["category dpmind"] = row["category"]


            # generate relationship id
            md5 = hashlib.md5(
                "".join(
                    [str(source_id), str(target_id)]
                ).encode("utf-8")
            ).hexdigest()

            # generate edge
            yield BioCypherEdge(
                relationship_id=md5,
                source_id=source_id,
                target_id=self._prefix_mrna(target_id),
                relationship_label="represses_dpmind",
                properties=properties,
            )
            
    def _set_types_and_fields(
        self,
        edge_types,
        edge_fields,
    ):

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in DpmindAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in DpmindAdapterMicrornaInhibitionEdgeField]

    @lru_cache(maxsize=None)
    
    def _prefix_mrna(self, string):
        return f"mrna:{string}"
    
    def _filter_input_dpmind(self):
        """Filter DPMIND using the filter_input_genome from GenomeAdapter

        Returns:
            dpmind_filtered: Filtered DataFrame of the original DPMIND input 
        """
    
        dpmind=pd.read_excel(self.file_path)
        dpmind['tar_OLN']=dpmind['tar_name'].str.split('.').str.get(0)
        dpmind['mir_name']=dpmind['mir_name'].replace('sly-miR171b','sly-miR171b-3p', regex=True)##Previously identify and should be rename
        
        # print('DPMIND interactions:', dpmind[['tar_OLN','mir_name']].drop_duplicates().shape[0])
        
        dpmind_filtered=GenomeAdapter(self.genome_path).filter_input_genome(dpmind,'tar_OLN')
        
        dpmind_filtered=MirbaseAdapter().filter_input_mirbase(dpmind_filtered,'mir_name')
        
        # print('DPMIND interactions after filtering:', dpmind_filtered[['tar_OLN','mir_name']].drop_duplicates().shape[0])
        
        return dpmind_filtered