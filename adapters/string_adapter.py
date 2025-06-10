from functools import lru_cache
import hashlib
from enum import Enum
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherEdge
from scripts.Uniprot_ID_mapping import *
from adapters.genome_adapter import GenomeAdapter
import os
from collections import defaultdict

logger.debug(f"Loading module {__name__}.")
    
class StringAdapterEdgeType(Enum):
    """
    Enum for the types of the String adapter
    """
    
    PROTEIN_PROTEIN_INTERACTION="interacts_with"
    
class StringAdapterMicrornaInhibitionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein protein interaction edges
    """
    
    SCORE="combined_score"
    
class StringAdapter:
    """
    STRING adapter. Import, filter and create links between proteins.

    Args:
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
        genome_path: str,
        combined_score: int,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self.genome_path = genome_path
        self.combined_score = combined_score
        self._set_types_and_fields(
            edge_types,
            edge_fields,
        )
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the given file and extract proteins.
        """
        logger.info("Preprocessing STRING data.")

        # load data
        self.data = self._filter_input_string()
            
    def get_edges(self):
        """

        Returns a generator of BioCypher edge objects (optionally
        BioCypherRelAsNode) for edge types specified in the adapter constructor.

        """

        logger.info("Generating edges.")

        # one row of the dataframe represents one edge
        for _, row in self.data.iterrows():
        
            # extract source and target
            source_id = row["OLN1"]
            target_id = row["OLN2"]
            #extract edge properties
            
            properties = {}
            
            if ( StringAdapterMicrornaInhibitionEdgeField.SCORE in self.edge_fields ):
                properties["score"] = row['combined_score']

            # generate relationship id
            md5 = hashlib.md5(
                "".join(
                    [str(source_id), str(target_id)]
                ).encode("utf-8")
            ).hexdigest()

            # generate edge
            yield BioCypherEdge(
                relationship_id=md5,
                source_id=self._prefix_protein(source_id),
                target_id=self._prefix_protein(target_id),
                relationship_label="interacts_with",
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
            self.edge_types = [type for type in StringAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in StringAdapterMicrornaInhibitionEdgeField]

    @lru_cache(maxsize=None)
    
    def _prefix_protein(self, string):
        return f"protein:{string}"
    
    def _create_uniprot_mapping(self):
        """ID conversion step: from UniProt ID to OLN ID. Filtering on the output OLN using the filter_input_genome from GenomeAdapter

        Returns:
            uniprot_mapping_filtered: DataFrame of the conversion with UniProt ID and OLN ID after filtering
        """
        
        logger.info("Mapping IDs.")
        
        string_aliases=pd.read_csv('download/string/4081.protein.aliases.v12.0.txt.gz.decomp', sep='\t')
        string_aliases['protein_id']=string_aliases['#string_protein_id'].str.split('.').str.get(1)
        string_protein_id=string_aliases['protein_id'].unique().tolist()
        job_id = submit_id_mapping(
            from_db="UniProtKB_AC-ID", to_db="Ensembl_Genomes", ids=string_protein_id
        )
            
        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results_map = get_id_mapping_results_search(link + "?format=tsv")
            
        uniprot_mapping=get_data_frame_from_tsv_results(results_map)
        uniprot_mapping['To_OLN']=uniprot_mapping['To'].str.split('.').str.get(0)
        
        # print('String proteins before filtering OLN:', uniprot_mapping['To_OLN'].unique().shape[0])
        
        uniprot_mapping_filtered=GenomeAdapter(self.genome_path).filter_input_genome(uniprot_mapping, 'To_OLN')
        
        # print('String proteins after filtering OLN:', uniprot_mapping_filtered['To_OLN'].unique().shape[0])
        
        return uniprot_mapping_filtered
    
    def _filter_input_string(self):
        """Filter the protein protein interaction file using the previous mapping and a threshold on combined score.

        Returns:
            string: Filtered DataFrame of the PPI input file with OLN and UniProt IDs
        """
        
        output_file='download/string/String_PPI_sly_filtered.txt'
        
        if not os.path.exists(output_file):
            logger.info("Creating file.")
            # Generate genome_ids directly from the function
            genome_ids = defaultdict(set)
            for _, row in self._create_uniprot_mapping().iterrows():
                from_id = row['From']
                to_id = row['To_OLN']
                genome_ids[from_id].add(to_id)

            # Filter the STRING database file and add genome IDs as additional columns
            seen_interactions = set()  # to keep track of unique interactions
            with open(output_file, 'w') as filtered_file:
                with open('download/string/4081.protein.links.v12.0.txt.gz.decomp', 'r') as string_file:
                    header = next(string_file).strip()  # Read and store the header
                    # Split the header into individual column names
                    column_names = header.split()
                    # Construct the new header with tabs as delimiters
                    header_with_genome_ids = '\t'.join(column_names) + '\tOLN1\tOLN2\n'
                    # Write the new header to the filtered file
                    filtered_file.write(header_with_genome_ids)

                    for line in string_file:
                        protein1, protein2, score = line.strip().split()
                        # Extract correct IDs from protein IDs in STRING file
                        protein1_id = protein1.split('.')[1] if '.' in protein1 else protein1
                        protein2_id = protein2.split('.')[1] if '.' in protein2 else protein2
                        # Check if both protein1 and protein2 are in the OLN dict. We keep the interaction only if both are present.
                        if protein1_id in genome_ids and protein2_id in genome_ids and int(score) > self.combined_score:  # combine score filter
                            interaction = tuple(sorted((protein1_id, protein2_id)))  # sort to ensure uniqueness
                            if interaction not in seen_interactions:
                                seen_interactions.add(interaction)
                                # Duplicate interactions for all combinations of genome IDs
                                for genome1_id in genome_ids[protein1_id]:
                                    for genome2_id in genome_ids[protein2_id]:
                                        # Adjust delimiter to '\t' in the output line
                                        line_with_genome_ids = f"{protein1}\t{protein2}\t{score}\t{genome1_id}\t{genome2_id}\n"
                                        filtered_file.write(line_with_genome_ids)
                
        else : 
            logger.info(f"File '{output_file}' already exists. Skipping creation.")
            
        string = pd.read_csv(output_file, sep='\t')
        
        # print('String PPI interaction:', string[['OLN1','OLN2']].drop_duplicates().shape[0])
        
        proteins = pd.concat([string['OLN1'], string['OLN2']])
        
        # print('String proteins after mapping and filtering:', proteins.unique().shape[0])
        
        return string