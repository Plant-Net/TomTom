from functools import lru_cache
import hashlib
from enum import Enum, auto
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherNode, BioCypherEdge
import subprocess
import os
import re

logger.debug(f"Loading module {__name__}.")

class MirbaseAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide
    """
    MICRORNA= auto()
    PRECURSOR = auto()
    
    
class MirbaseAdapterMicrornaField(Enum):
    """
    Define possible fields the adapter can provide for micrornas
    """
    
    MICRORNA_ID="Mature ID"
    ACCESSION="Mature Accession"
    SEQUENCE="Mature Sequence"

class MirbaseAdapterPrecursorField(Enum):
    """
    Define possible fields the adapter can provide for precursors
    """
    
    PRECURSOR_ID="Precursor ID"
    ACCESSION="Precursor Accession"
    SEQUENCE="Precursor Sequence"
    
class MirbaseAdapterEdgeType(Enum):
    """
    Enum for the edge types of the microrna adapter
    """
    
    MICRORNA_PRODUCTION="matures_to"
    
class MirbaseAdapterMicrornaProductionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for microrna production edges
    """
    
    START_POSITION = "Start Position"
    END_POSITION = "End Position"
    
class MirbaseAdapter:
    """
    MiRBase adapter. Import and create microRNA and precursor nodes and edges for creating a
    knowledge graph backbone.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
    """

    def __init__(
        self,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self._set_types_and_fields(
            node_types,
            node_fields,
            edge_types,
            edge_fields,
        )
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the given file and extract micrornas and precursors.
        """
        
        self.data = self.read_parse_miRBase()
        
        # extract precursors (unique entities of `precursor` column)
        self.precursors = self.data[["Precursor ID","Precursor Accession","Precursor Sequence"]].drop_duplicates()

        # extract micrornas (unique entities of `microrna` column)
        self.micrornas = self.data[["Mature ID","Mature Accession","Mature Sequence"]].drop_duplicates()
        
    def get_nodes(self):
        """

        Returns a generator of BioCypher node objects for node types specified
        in the adapter constructor.

        Returns:
            Generator of BioCypher node objects.

        """

        logger.info("Generating nodes.")

        for _, row in self.precursors.iterrows():
            node_id=row["Precursor ID"]
            accession=row["Precursor Accession"]
            sequence=row["Precursor Sequence"]
            properties = {
                "name":node_id,
                "accession":accession,
                "sequence":sequence
            }
        
            yield BioCypherNode(
                node_id=node_id,
                node_label="precursor",
                properties=properties,
            )


        for _, row in self.micrornas.iterrows():
            node_id=row["Mature ID"]
            accession=row["Mature Accession"]
            sequence=row["Mature Sequence"]
            properties = {
                "name":node_id,
                "accession":accession,
                "sequence":sequence
            }
        
            yield BioCypherNode(
                node_id=node_id,
                node_label="microRNA",
                preferred_id="microrna_id",
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
            source_id = row["Precursor ID"]
            target_id = row["Mature ID"]
            
            #extract edge properties
            
            properties = {}
            
            if ( MirbaseAdapterMicrornaProductionEdgeField.START_POSITION in self.edge_fields ):
                properties["start position"] = row["Start Position"]
                
            if ( MirbaseAdapterMicrornaProductionEdgeField.END_POSITION in self.edge_fields ):
                properties["end position"] = row["End Position"]



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
                target_id=target_id,
                relationship_label="matures_to",
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
            self.node_types = [type for type in MirbaseAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in MirbaseAdapterMicrornaField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in MirbaseAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in MirbaseAdapterMicrornaProductionEdgeField]

    @lru_cache(maxsize=None)
    
    def extract_sol_mirnas(self):
        with open('download/mirbase/miRNA.dat', 'r') as f:
            data = f.read()

        # Extract blocks corresponding to each miRNA
        miRNAs = re.findall(r'^\s*ID.*?\n//', data, re.S | re.M)

        sol_mirnas = {}

        for mirna in miRNAs:
            # Ensure it belongs to Solanum lycopersicum
            if "Solanum lycopersicum" in mirna:
                # Extract the precursor miRNA ID
                precursor_match = re.search(r'ID\s+(\S+)', mirna)
                if precursor_match:
                    precursor_id = precursor_match.group(1)
                    
                    # Extract the precursor accession number
                    accession_match = re.search(r'AC\s+(\S+);', mirna)
                    precursor_accession = accession_match.group(1) if accession_match else None
                    
                    # Extract the mature miRNAs blocks within this precursor entry
                    mature_mirnas = []
                    for match in re.finditer(r'FT\s+miRNA\s+(\d+..\d+)\nFT\s+/accession="(\S+)"\nFT\s+/product="([^"]+)"', mirna):
                        position_range = match.group(1)
                        start, end = map(int, position_range.split('..'))
                        mature_accession = match.group(2)
                        product = match.group(3)
                        mature_mirnas.append({
                            "accession": mature_accession,
                            "product": product,
                            "start": start,
                            "end": end
                        })

                    # Store precursor and its corresponding mature miRNAs in the final dictionary
                    if mature_mirnas:
                        sol_mirnas[precursor_id] = {
                            'accession': precursor_accession,
                            'mature_products': mature_mirnas
                        }

        return sol_mirnas


    def extract_mature_sequences(self, sol_mirnas):
        mature_sequences = {}

        with open('download/mirbase/mature.fa', 'r') as f:
            lines = f.readlines()

        current_mature_id = None
        current_sequence = []

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # Parse the header to extract the mature miRNA ID (product name)
                header = line.split()
                mature_id = header[0][1:]  # Remove '>'
                
                # Check if the mature ID is in the sol_mirnas
                found = False
                for precursor_id, info in sol_mirnas.items():
                    for mature in info['mature_products']:
                        if mature['product'] == mature_id:
                            current_mature_id = mature['product']
                            found = True
                            break
                    if found:
                        break

                # If we had a previous mature ID, store its sequence
                if current_mature_id and current_sequence:
                    mature_sequences[current_mature_id] = ''.join(current_sequence)

                # Reset for the next mature miRNA
                current_sequence = []
            else:
                # Accumulate the sequence lines for the current mature miRNA
                current_sequence.append(line)

        # Store the last sequence after the loop
        if current_mature_id and current_sequence:
            mature_sequences[current_mature_id] = ''.join(current_sequence)

        return mature_sequences


    def extract_precursor_sequences(self, sol_mirnas):
        precursor_sequences = {}

        with open('download/mirbase/hairpin.fa', 'r') as f:
            lines = f.readlines()

        current_precursor_id = None
        current_sequence = []

        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # Parse the header to extract the precursor miRNA ID
                header = line.split()
                precursor_id = header[0][1:]  # Remove '>'
                
                # Check if the precursor ID exists in the sol_mirnas dictionary
                if precursor_id in sol_mirnas:
                    current_precursor_id = precursor_id

                # If we had a previous precursor ID, store its sequence
                if current_precursor_id and current_sequence:
                    precursor_sequences[current_precursor_id] = ''.join(current_sequence)

                # Reset for the next precursor miRNA
                current_sequence = []
            else:
                # Accumulate the sequence lines for the current precursor miRNA
                current_sequence.append(line)

        # Store the last sequence after the loop
        if current_precursor_id and current_sequence:
            precursor_sequences[current_precursor_id] = ''.join(current_sequence)

        return precursor_sequences



    def sol_mirnas_to_dataframe(self, sol_mirnas, precursor_sequences, mature_sequences):
        """
        Convert sol_mirnas dictionary to a DataFrame, including precursor and mature miRNA sequences.
        
        Parameters:
        - sol_mirnas: Dictionary containing miRNA information.
        - precursor_sequences: Dictionary mapping precursor IDs to their sequences.
        - mature_sequences: Dictionary mapping mature product IDs to their sequences.
        
        Returns:
        - pd.DataFrame: A DataFrame containing all relevant information, including sequences.
        """
        
        # Create a list to store rows for DataFrame
        rows = []
        
        for precursor_id, info in sol_mirnas.items():
            # Extract precursor information
            precursor_accession = info['accession']
            precursor_seq = precursor_sequences.get(precursor_id, None)  # Get precursor sequence
            
            # Extract mature miRNAs
            for mature in info['mature_products']:
                mature_seq = mature_sequences.get(mature['product'], None)  # Get mature sequence
                
                row = {
                    'Precursor ID': precursor_id,
                    'Precursor Accession': precursor_accession,
                    'Precursor Sequence': precursor_seq,
                    'Mature ID': mature['product'],
                    'Mature Accession': mature['accession'],
                    'Mature Sequence': mature_seq,
                    'Start Position': mature['start'],
                    'End Position': mature['end']
                }
                rows.append(row)
        
        # Create DataFrame
        mirbase_df = pd.DataFrame(rows)
        return mirbase_df
    
    def read_parse_miRBase(self):
        sol_mirnas = self.extract_sol_mirnas()
        
        precursor_sequences = self.extract_precursor_sequences(sol_mirnas)
        
        mature_sequences = self.extract_mature_sequences(sol_mirnas)
        
        miRBase_df = self.sol_mirnas_to_dataframe(sol_mirnas, precursor_sequences, mature_sequences)
        
        return miRBase_df
    
    def filter_input_mirbase(self,df_to_filter,column):
        """Filter microRNAs of the input database 

        Args:
            df_to_filter (pandas DataFrame): DataFrame to filter
            column (str): df_to_filter column of the dataframe on which to perform the filter

        Returns:
            Filtered_df (pandas DataFrame): Filtered dataframe of the input
        """
        
        logger.info("Filtering...")
    
        miRBase=self.read_parse_miRBase()
        set_mirbase=set(miRBase['Mature ID'].unique())
        set_df=set(df_to_filter[column].unique())
        
        to_filter_out=set_df.difference(set_mirbase)
        
        # print("microRNA lost:", len(to_filter_out))
        
        Filtered_df=df_to_filter[~df_to_filter[column].isin(to_filter_out)]
        
        return Filtered_df