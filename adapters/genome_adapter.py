from functools import lru_cache
import hashlib
from enum import Enum, auto
from typing import Optional
import pandas as pd
from biocypher._logger import logger
from biocypher._create import BioCypherNode, BioCypherEdge
import re

logger.debug(f"Loading module {__name__}.")

class GenomeAdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide
    """
    GENE = auto()
    TRANSCRIPTION_FACTOR = auto()
    MRNA = auto()
    PROTEIN = auto()
    

class GenomeAdapterGeneField(Enum):
    """
    Define possible fields the adapter can provide for genes
    """
    
    GENE_SYMBOL="symbol"
    DESCRIPTION = "Description"
    ALIAS = "Alias"
    
class GenomeAdapterTfField(Enum):
    """
    Define possible fields the adapter can provide for transcription factors
    """
    
    GENE_SYMBOL="symbol"
    DESCRIPTION = "Description"
    FAMILY= "Family"
    ALIAS = "Alias"
    
class GenomeAdapterMrnaField(Enum):
    """
    Define possible fields the adapter can provide for mrnas
    """
    
    MRNA_SYMBOL="symbol"
    DESCRIPTION = "Description"
    ALIAS = "Alias"
    
class GenomeAdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins
    """
    
    PROTEIN_SYMBOL="symbol"
    DESCRIPTION = "Description"
    ALIAS = "Alias"
        
class GenomeAdapterEdgeType(Enum):
    """
    Enum for the edge types of the genome adapter
    """
    
    MRNA_PRODUCTION="transcripted_to"
    PROTEIN_PRODUCTION="translated_to"
    
class GenomeAdapterMrnaProductionEdgeField(Enum):
    """
    Define possible fields the adapter can provide for mrna production edges
    """
    
    CHR = "chr"
    START = "start"
    END = "end"
    STRAND = "strand"
    

class GenomeAdapter:
    """
    Genome adapter. Import the genome and write nodes and edges for creating a
    knowledge graph backbone.

    Args:
        node_types: List of node types to include in the result.
        node_fields: List of node fields to include in the result.
        edge_types: List of edge types to include in the result.
        edge_fields: List of edge fields to include in the result.
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
        self._genome_parsed = self._parse_gff_genome()
        self._preprocess_data()
    
        
    def _preprocess_data(self):
        """
        Load the data from the given GFF and extract genes, mrnas and proteins.
        """
        logger.info("Preprocessing Genome data.")

        # load data
        self.data = self._create_genome_info()

        # extract genes (unique entities of `target` column)
        self.genes = self.data[["OLN","Description", "Family", "Alias"]]

        # extract mrnas (unique entities of `source` column)
        self.mrnas = self.data[["OLN","Description", "Alias"]]
        
        #extract proteins
        self.proteins = self.data[["OLN","Description", "Alias"]]
        
    def get_nodes(self):
        """

        Returns a generator of BioCypher node objects for node types specified
        in the adapter constructor.

        Returns:
            Generator of BioCypher node objects.

        """

        logger.info("Generating nodes.")
     
        for _, row in self.genes.iterrows():
            node_id = row["OLN"]
            description = row["Description"]
            alias = row["Alias"]
            properties = {
                "name": node_id,
                "description": description,
                "alias": alias
            }
            
            if pd.notna(row['Family']):
                family = row['Family']
                properties["family"] = family
                yield BioCypherNode(
                    node_id=self._prefix_gene(node_id),
                    node_label="transcription factor",
                    properties=properties,
                )
            else:
                yield BioCypherNode(
                    node_id=self._prefix_gene(node_id),
                    node_label="gene",
                    properties=properties,
                )
            
        for _, row in self.mrnas.iterrows():
            node_id=row["OLN"]
            description=row["Description"]
            alias=row["Alias"]
            properties = {
                "name":node_id,
                "description": description,
                "alias": alias
            }
    
            yield BioCypherNode(
                node_id=self._prefix_mrna(node_id),
                node_label="mrna",
                properties=properties,
            )
            
        for _, row in self.proteins.iterrows():
            node_id=row["OLN"]
            description=row["Description"]
            alias=row["Alias"]
            properties = {
                "name":node_id,
                "description": description,
                "alias": alias
            }
    
            yield BioCypherNode(
                node_id=self._prefix_protein(node_id),
                node_label="protein",
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
            target_id = row["OLN"]
            
            #extract edge properties
            
            properties = {}
            
            if ( GenomeAdapterMrnaProductionEdgeField.CHR in self.edge_fields ):
                properties["chr"] = row["chr"]
            
            if ( GenomeAdapterMrnaProductionEdgeField.START in self.edge_fields ):
                properties["start"] = row["start"]
                
            if ( GenomeAdapterMrnaProductionEdgeField.END in self.edge_fields ):
                properties["end"] = row["end"]
            
            if ( GenomeAdapterMrnaProductionEdgeField.STRAND in self.edge_fields ):
                properties["strand"] = row["strand"]


            # generate relationship id
            md5 = hashlib.md5(
                "".join(
                    [str(source_id), str(target_id)]
                ).encode("utf-8")
            ).hexdigest()

            # generate edge
    
            # gene to mrna
            yield BioCypherEdge(
                relationship_id=md5,
                source_id=self._prefix_gene(source_id),
                target_id=self._prefix_mrna(target_id),
                relationship_label="transcripted_to",
                properties=properties,
            )
            
            # mrna to protein
            yield BioCypherEdge(
                relationship_id=md5,
                source_id=self._prefix_mrna(source_id),
                target_id=self._prefix_protein(target_id),
                relationship_label="translated_to",
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
            self.node_types = [type for type in GenomeAdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [field for field in GenomeAdapterMrnaField]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in GenomeAdapterEdgeType]
            
        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in GenomeAdapterMrnaProductionEdgeField]

    @lru_cache(maxsize=None)
    #Define prefix to distinguish gene, mrna and protein
    def _prefix_gene(self, string):
        return f"gene:{string}"
    def _prefix_mrna(self, string):
        return f"mrna:{string}"
    def _prefix_protein(self, string):
        return f"protein:{string}"
    
    def _parse_gff_genome(self):
        """Read and parse one GFF file to get an edge list of the genome

        Returns:
            Genome_list (pandas DataFrame): DataFrame of the genome in an edge list format with gene, mrna, chr, start, end and strand
        """
        logger.info("Parsing genome GFF file.")
        gene_mrna_pairs = []

        with open(self.genome_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                attributes = fields[8]
                attributes_dict = dict(attribute.split('=') for attribute in attributes.split(';') if '=' in attribute)

                if feature_type == 'gene':
                    gene_id = attributes_dict.get('ID', '')
                    gene_id = gene_id.split(':')[1] if ':' in gene_id else gene_id
                elif feature_type == 'mRNA':
                    mrna_id = attributes_dict.get('ID', '')
                    mrna_id = mrna_id.split(':')[1] if ':' in mrna_id else mrna_id
                    parent_gene_id = attributes_dict.get('Parent', '')
                    parent_gene_id = parent_gene_id.split(':')[1] if ':' in parent_gene_id else parent_gene_id
                    gene_mrna_pairs.append((parent_gene_id, mrna_id, fields[0], fields[3], fields[4], fields[6]))
                    
        Genome_list=pd.DataFrame(gene_mrna_pairs, columns=['gene', 'mrna', 'chr', 'start', 'end','strand'])
        Genome_list['OLN'] = Genome_list['gene'].str.split('.').str.get(0)
        
        return Genome_list
   
    def _read_parse_description(self):
        """Read the genome description file
        
        Returns: 
            description_df (pandas DataFrame): dataframe containing gene ID and associated description
        """
            # Read the file and split lines
        with open('download/genome/ITAG4.1_descriptions.txt', "r") as file:
            lines = file.readlines()

        # Define a regular expression pattern to extract the two columns
        pattern = r'^(.*?)\s(.*?)\s*$'

        # Initialize lists to store data
        ids = []
        descriptions = []

        # Iterate through each line and extract data
        for line in lines:
            match = re.match(pattern, line)
            if match:
                ids.append(match.group(1))
                descriptions.append(match.group(2))

        # Create pandas DataFrame
        description_df = pd.DataFrame({
            'ID': ids,
            'Description': descriptions
        })
        
        description_df['OLN'] = description_df['ID'].str.split('.').str.get(0)
        description_df['Description'] = description_df['Description'].str.replace("'","''")
        
        return description_df
    
    def _merge_info(self):
        """Merge the genome informations (genes, position,...) with description

        Returns:
            info_genome (pandas DataFrame): Merge dataframe
        """
        genome = self._genome_parsed
        description = self._read_parse_description()
        
        info_genome = pd.merge(genome, description, on='OLN')
        
        return info_genome
    
    def filter_input_genome(self, df_to_filter, column):
        """Filter an input dataframe coming from other databases

        Args:
            df_to_filter (pandas DataFrame): Input dataframe to filter. Must contain an OLN column to identify entity
            column (str): df_to_filter column on which to perform filtering

        Returns:
            Filtered_df (pandas DataFrame): Output filtered dataframe
        """
        logger.info("Filtering...")
        
        Genome=self._genome_parsed
        set_genome=set(Genome['OLN'].unique())
        set_df=set(df_to_filter[column].unique())
        to_filter_out=set_df.difference(set_genome)
        
        # logger.info('Genomic entities lost:', len(to_filter_out))
        # logger.info('First 5 genomic entities lost:', list(to_filter_out)[:5])
        
        
        Filtered_df=df_to_filter[~df_to_filter[column].isin(to_filter_out)]
        
        return Filtered_df
    
    def _filter_input_planttfdb(self):
        """Filter the planttfdb input. Allow the idenficiation of transcription factors

        Returns:
            planttfdb_filtered (pandas DataFrame): DataFrame of the database after filtering containing gene ID and family of transcription factors
        """
        
        planttfdb=pd.read_csv('download/planttfdb/Sly_TF_list.txt.gz.decomp', sep='\t')
        planttfdb['OLN']=planttfdb['TF_ID'].str.split('.').str.get(0)
        
        # logger.info('Total interaction:', planttfdb.shape[0])
        # logger.info('Total unique TF:', planttfdb['OLN'].unique().shape[0])
        
        planttfdb_filtered=self.filter_input_genome(planttfdb, 'OLN')
        planttfdb_filtered=planttfdb_filtered[['OLN','Family']]
        
        # logger.info('Total interaction after filtering:', planttfdb_filtered.shape[0])
        # logger.info('Total unique TF after filtering:', planttfdb_filtered['OLN'].unique().shape[0])
        
        return planttfdb_filtered
    
    def _read_format_synonyms(self):
        """Read and format the synonyms file to get a dataframe with OLN and associated synonyms

        Returns:
            synonyms_df (pandas DataFrame): DataFrame containing OLN and associated synonyms
        """
        synonyms_df = pd.read_csv('download/biomart/biomart_slycopersicum_synonyms.tsv', sep='\t', header=None)
        synonyms_df = synonyms_df.rename(columns={0: 'OLN'})
        synonyms_df = synonyms_df[synonyms_df["OLN"].str.startswith('Solyc')]
        synonyms_df['OLN'] = synonyms_df['OLN'].str.split('.').str.get(0)
        
        long_synonyms_df = (synonyms_df.melt(id_vars='OLN', value_vars=[1,2], value_name='Alias').dropna().drop_duplicates(["OLN","Alias"]))
        formatted_synonyms_df = (long_synonyms_df.groupby("OLN")["Alias"].agg(list).reset_index())
        
        return formatted_synonyms_df
    
    def _create_genome_info(self):
        """Create the genome information from previous dataframes by merging the information (genes, positions, description, ...) with planttfdb and aliases

        Returns:
            genome_info_tf (pandas DataFrame): DataFrame of the genome with all the information needed to construct the backbone of the graph.
        """
        genome_info = self._merge_info()
        planttfdb = self._filter_input_planttfdb()
        Aliases = self._read_format_synonyms()
        
        genome_info_tf = pd.merge(planttfdb, genome_info, on='OLN', how='right')
        genome_info_tf_alias = pd.merge(genome_info_tf, Aliases, on='OLN', how='left')
        genome_info_tf_alias['Alias'] = genome_info_tf_alias['Alias'].fillna('NA')
        
        return genome_info_tf_alias.drop_duplicates(subset=['OLN','Description','Family'])
        