# 🍅 Tomato Knowledge Graph 🍅
*Solanum lycopersicum* (tomato) is a plant of major agronomic interest and an increasingly studied organism.  
However, knowledge on the tomato is widely spread across different databases. Bringing together information on this organism in one place could help many biologists to speed up their understanding.  
Here I developed a knowledge graph (KG) on *S. lycopersicum* species using [BioCypher](https://github.com/biocypher/biocypher).
- [🍅 Tomato Knowledge Graph 🍅](#-tomato-knowledge-graph-)
  - [Input Databases 📚](#input-databases-)
  - [Installation ⚙️](#installation-️)
  - [Docker 🐳](#docker-)
  - [Reference :writing\_hand:](#reference-writing_hand)
## Input Databases 📚
The KG is composed of several input databases as described in the following table : 

| Database      | Description |
| :---:        |    :----:   |
| [Sol Genomics Network](https://solgenomics.net/)      | Genome       |
| [miRBase](https://mirbase.org/)   |  microRNA and precursor       |
| [PlantTFDB](https://planttfdb.gao-lab.org/)   |  TF identification       |
| [PlantRegMap](https://plantregmap.gao-lab.org/)   |  TF-target interaction       |
| [TarDB](http://www.biosequencing.cn/TarDB/)   |  microRNA-transcript interaction       |
| [DPMIND](https://cbi.njau.edu.cn/DPMIND/)   |  microRNA-transcript interaction       |
| [PNRD](https://structuralbiology.cau.edu.cn/PNRD/index.php)   |  microRNA-transcript interaction        |
| [STRING](https://string-db.org/)  |  protein-protein interaction      |
| [Planteome](https://planteome.org/)  |  term associated to gene       |
| [Mercator4](https://www.plabipd.de/mercator_main.html)   |  pathway associated to gene       |
| [OMA](https://omabrowser.org/oma/home/) | Gene - *A.thaliana* gene  |

## Installation ⚙️
Once you clone the repository, you can install the dependencies using poetry:

```bash
poetry install
```
Then, you should be able to create the knowledge graph by first downloading all the databases. 
The databases must be downloaded before creating the graph.

```bash
poetry shell
python scripts/download_databases.py
python create_knowledge_graph.py
```
If everything runs smoothly, you can run the Docker 🐳

NB: You can exit the poetry shell just by typing ```exit```

## Docker 🐳

After downloading the files and make sure the graph can be built, we can start the Neo4j database with the docker:

```bash
docker compose up -d
```
You can connect and browse the Neo4j instance at localhost:7474. No authentification is needed, just press connect.

To shutdown the docker : 
```bash
docker compose down -v
```

## Reference :writing_hand: 
You can find all the details on the graph construction in the [publication](https://doi.org/10.1111/nph.71382)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15630284.svg)](https://doi.org/10.5281/zenodo.15630284)
