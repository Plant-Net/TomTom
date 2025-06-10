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
    
class PnrdAdapterEdgeType(Enum):
    """
    Enum for the types of the MTI adapter
    """
    
    MICRORNA_INHIBITION="represses_pnrd"
    
class PnrdAdapterMicrornaInhibitionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for microrna inhibition edges
    """
    
    EXPECTION_PNRD="Expection"
    METHOD_PNRD="method"
    
class PnrdAdapter:
    """
    PNRD adapter. Import, filter and create links between microRNAs and targets.

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
        logger.info("Preprocessing PNRD data.")

        # load data
        self.data = self._filter_input_pnrd()
            
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.data.iterrows():
            # extract source and target
            source_id = row["ncRNA name"]
            target_id = row["target OLN"]
            
            #extract edge properties
            
            properties = {}
            
            if ( PnrdAdapterMicrornaInhibitionEdgeField.EXPECTION_PNRD in self.edge_fields ):
                properties["expection pnrd"] = row["Expection"]
            
            if ( PnrdAdapterMicrornaInhibitionEdgeField.METHOD_PNRD in self.edge_fields ):
                properties["method pnrd"] = row["method"]
            

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
                relationship_label="represses_pnrd",
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
            self.edge_types = [type for type in PnrdAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in PnrdAdapterMicrornaInhibitionEdgeField]

    @lru_cache(maxsize=None)
    
    def _prefix_mrna(self, string):
        return f"mrna:{string}"
    
    def _filter_input_pnrd(self):
        """Filter PNRD using the filter_input_genome from GenomeAdapter

        Returns:
            pnrd_filtered_final: Filtered DataFrame of the original PNRD input 
        """
    
        pnrd=pd.read_csv(self.file_path, sep='\t', header=None)
        ##Create header based on the web page
        pnrd.rename(columns={0:'ncRNA name',1:'target name',2:'Expection',3:'ncRNA Alignment start',4:'ncRNA Aligment end',5:'target Alignment start',6:'target Alignment end',7:'ncRNA sequence',8:'target sequence',9:'target description',10:'method',11:'link',12:'target type'}, inplace=True)
        pnrd['target OLN']=pnrd['target name'].str.split('.').str.get(0)
        pnrd['ncRNA sequence']= pnrd['ncRNA sequence'].apply(lambda x: x[::-1]) #Need to invert the sequence to match miRBase info
        
        # print('PNRD interactions:', pnrd[['target OLN','ncRNA name']].drop_duplicates().shape[0])
        
        pnrd_filtered=GenomeAdapter(self.genome_path).filter_input_genome(pnrd,'target OLN')
        
        ##PNRD had some problem with miRBase ID. Some of the id reflect precursor and not mature microRNA. By applying a merge on the sequence column, we can correct the issue.
        
        miRBase=MirbaseAdapter().read_parse_miRBase()
        
        set_pnrd_mirna=set(pnrd_filtered['ncRNA name'].unique())
        set_mirbase=set(miRBase['Mature ID'].unique())
        
        problem_mirna_pnrd=pnrd_filtered[pnrd_filtered['ncRNA name'].isin(set_pnrd_mirna.difference(set_mirbase))]
        
        match_seq=pd.merge( miRBase,problem_mirna_pnrd,left_on='Mature Sequence', right_on='ncRNA sequence', how='inner' ) #Perform a merge on sequence to have 100% identity
        
        dict_replacement=dict(match_seq[['ncRNA name','Mature ID']].drop_duplicates().values)
        
        pnrd_filtered.loc[:,'ncRNA name']=pnrd_filtered['ncRNA name'].replace(dict_replacement)
        
        pnrd_filtered_final=MirbaseAdapter().filter_input_mirbase(pnrd_filtered, 'ncRNA name')
        
        # print('PNRD interactions after filtering:', pnrd_filtered_final[['target OLN','ncRNA name']].drop_duplicates().shape[0])
        
        return pnrd_filtered_final
