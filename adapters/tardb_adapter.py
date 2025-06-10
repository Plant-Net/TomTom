from functools import lru_cache
import hashlib
from enum import Enum
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherEdge
import os
import subprocess
from adapters.genome_adapter import GenomeAdapter
from adapters.miRBase_adapter import MirbaseAdapter

logger.debug(f"Loading module {__name__}.")
    
class TardbAdapterEdgeType(Enum):
    """
    Enum for the types of the MTI adapter
    """
    
    MICRORNA_INHIBITION="represses_tardb"
    
class TardbAdapterMicrornaInhibitionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for microrna inhibition edges
    """
    
    SCORE_TARDB = "Score"
    TOTAL_MISPAIR_TARDB = "Total_Mispair"
    SEED_MISPAIR_TARDB = "Seed_Mispair"
    PREDICTED_CLEAVAGE_TARDB = "Pred_Cleav"
    
class TardbAdapter:
    """
    TarDB adapter. Import, filter and create links between microRNAs and targets.

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
        logger.info("Preprocessing TarDB data.")

        # load data
        self.data = self._filter_input_tardb()
            
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.data.iterrows():
            # extract source and target
            source_id = row["miRNA_ID"]
            target_id = row["Target_OLN"]
            
            #extract edge properties
            
            properties = {}
            
            if ( TardbAdapterMicrornaInhibitionEdgeField.SCORE_TARDB in self.edge_fields ):
                properties["score tardb"] = row["Score"]
            
            if ( TardbAdapterMicrornaInhibitionEdgeField.TOTAL_MISPAIR_TARDB in self.edge_fields ):
                properties["total mispair tardb"] = row["Total_Mispair"]
                
            if ( TardbAdapterMicrornaInhibitionEdgeField.SEED_MISPAIR_TARDB in self.edge_fields ):
                properties["seed mispair tardb"] = row["Seed_Mispair"]
            
            if ( TardbAdapterMicrornaInhibitionEdgeField.PREDICTED_CLEAVAGE_TARDB in self.edge_fields ):
                properties["predicted cleavage tardb"] = row["Pred_Cleav"]


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
                relationship_label="represses_tardb",
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
            self.edge_types = [type for type in TardbAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in TardbAdapterMicrornaInhibitionEdgeField]

    @lru_cache(maxsize=None)
    
    def _prefix_mrna(self, string):
        return f"mrna:{string}"
    
    def _filter_input_tardb(self):
        """Filter TarDB using the filter_input_genome from GenomeAdapter

        Returns:
            tardb_filtered (pandas DataFrame): Filtered DataFrame of the original TarDB input 
        """
        
        #File path
        tardb_file = os.path.join(self.file_path,"sly","sly.cons")
        # Define the output file path
        output_file_path = os.path.join(self.file_path, "sly", "TarDB_MTI_sly.txt")
        
        # Define the bash script to extract data from the file
        bash_script = f"""
        awk 'NR == 1' {tardb_file}
        awk 'NR>1' {tardb_file} | awk 'NR % 5 == 2'
        """
        # Execute the bash script
        output = subprocess.run(bash_script, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Check for errors
        if output.returncode == 0:
        # Write the output to the file in the same directory
            with open(output_file_path, "w") as f:
                f.write(output.stdout)
        else:
            logger.error("Error occurred:", output.stderr)
        # Read the output
        tardb=pd.read_csv(output_file_path, sep='\t')
        tardb['Target_OLN']=tardb['Target_ID'].str.split('.').str.get(0)
        
        # print('Tardb interactions:', tardb[['Target_OLN','miRNA_ID']].drop_duplicates().shape[0])
        
        tardb_filtered=GenomeAdapter(self.genome_path).filter_input_genome(tardb,'Target_OLN')
        
        tardb_filtered=MirbaseAdapter().filter_input_mirbase(tardb_filtered,'miRNA_ID')
        
        # print('Tardb filtered interactions:', tardb_filtered[['Target_OLN','miRNA_ID']].drop_duplicates().shape[0])
        
        return tardb_filtered