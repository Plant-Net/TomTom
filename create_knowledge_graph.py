from biocypher import BioCypher
from adapters.genome_adapter import GenomeAdapter
from adapters.tardb_adapter import TardbAdapter
from adapters.miRBase_adapter import MirbaseAdapter
from adapters.dpmind_adapter import DpmindAdapter
from adapters.plantregmap_adapter import PlantregmapAdapter
from adapters.pnrd_adapter import PnrdAdapter
from adapters.string_adapter import StringAdapter
from adapters.planteome_adapter import PlanteomeAdapter
from adapters.kegg_adapter import KeggAdapter
from adapters.oma_adapter import OmaAdapter

bc=BioCypher()

genome = 'download/genome/ITAG4.1_gene_models.gff'
plantregmap= 'download/plantregmap/regulation_merged_Sly.txt'
dpmind = 'download/dpmind/Solanum_lycopersicum.xlsx'
tardb = 'download/tardb/sly.zip.unzip'
pnrd = 'database_downloads/database_downloads/PNRD/PNRD_sly_targets.txt'

#Building the genome backbone
genome_adapter = GenomeAdapter(genome_path = genome)

#Building the microRNA backbone
mirbase_adapter = MirbaseAdapter()

#Adding the TF regulation
TF_adapter=PlantregmapAdapter(plantregmap,genome_path=genome)

#Adding MTI
dpmind_adapter= DpmindAdapter(dpmind, genome_path=genome)

tardb_adapter = TardbAdapter(tardb, genome_path=genome)

pnrd_adapter = PnrdAdapter(pnrd, genome_path=genome)

# Adding String
string_adapter = StringAdapter(genome_path=genome, combined_score=800)

#Adding planteome terms
planteome_adapter = PlanteomeAdapter(genome_path=genome)

#Adding KEGG Pathways
kegg_adapter = KeggAdapter(genome_path=genome)

#Adding OMA
oma_adapter = OmaAdapter(genome_path=genome)    

#Adding the nodes and edges to the knowledge graph
adapters_for_nodes = [genome_adapter, mirbase_adapter, planteome_adapter, kegg_adapter, oma_adapter]
adapters_for_edges = [genome_adapter, mirbase_adapter, TF_adapter, dpmind_adapter, tardb_adapter, 
                      pnrd_adapter, string_adapter, planteome_adapter, kegg_adapter, oma_adapter]

nodes = (node for adapter in adapters_for_nodes for node in adapter.get_nodes())
edges = (edge for adapter in adapters_for_edges for edge in adapter.get_edges())

bc.write_nodes(nodes)
bc.write_edges(edges)

bc.write_import_call()

# bc.write_schema_info(as_node=True)