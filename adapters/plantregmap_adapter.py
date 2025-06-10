from functools import lru_cache
import hashlib
from enum import Enum
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherEdge
from adapters.genome_adapter import GenomeAdapter

logger.debug(f"Loading module {__name__}.")
    
class PlantregmapAdapterEdgeType(Enum):
    """
    Enum for the edge types of the plantregmap adapter
    """
    
    TF_REGULATION="regulates"
    
class PlantregmapAdapterMicrornaInhibitionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for TF Regulation edges
    """
    
    EVIDENCE="evidence"
    
class PlantregmapAdapter:
    """
    PlantRegMap adapter. Import, filter and create links between transcription factors and targets.

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
        Load the data from the given file and filter the entity.
        """
        logger.info("Preprocessing PlantRegMap data.")

        # load data
        self.data = self._filter_input_plantregmap()
            
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        
        # one row of the dataframe represents one edge
        
        for _, row in self.data.iterrows():
            # extract source and target
            source_id = row["Source_OLN"]
            target_id = row["Target_OLN"]
            
            #extract edge properties
            
            properties = {}
            
            if ( PlantregmapAdapterMicrornaInhibitionEdgeField.EVIDENCE in self.edge_fields ):
                properties["evidence"] = row["Evidence"]
            

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
                target_id=self._prefix_gene(target_id),
                relationship_label="regulates",
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
            self.edge_types = [type for type in PlantregmapAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in PlantregmapAdapterMicrornaInhibitionEdgeField]

    @lru_cache(maxsize=None)
    
    def _prefix_gene(self, string):
        return f"gene:{string}"
    
    def _filter_input_plantregmap(self):
        """Filter PlantRegMap using the filter_input_genome from GenomeAdapter

        Returns:
            plantregmap_filtered (pandas DataFrame): Filtered DataFrame of the original PlantRegMap input 
        """
        
        logger.info("Filtering PlantRegMap.")
    
        plantregmap=pd.read_csv(self.file_path, sep='\t', header=None)
        plantregmap['Source_OLN']=plantregmap[0].str.split('.').str.get(0)
        plantregmap['Target_OLN']=plantregmap[2].str.split('.').str.get(0)
        plantregmap.rename(columns={4:'Evidence'},inplace=True)
        
        # print('PlantRegMap interactions:', plantregmap[['Source_OLN','Target_OLN']].drop_duplicates().shape[0])
        
        plantregmap_filtered=GenomeAdapter(self.genome_path).filter_input_genome(plantregmap, 'Source_OLN')
        
        plantregmap_filtered=GenomeAdapter(self.genome_path).filter_input_genome(plantregmap_filtered, 'Target_OLN')
        
        # print('PlantRegMap interactions after filtering:', plantregmap_filtered[['Source_OLN','Target_OLN']].drop_duplicates().shape[0])
        
        return plantregmap_filtered