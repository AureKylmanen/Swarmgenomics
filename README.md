# SwarmGenomics


Advances in sequencing technologies have led to a boom in genomic projects, making genomic data readily available. This Citizen Science Project aims to demonstrate how a single genome can be analysed with simple computational methods to provide valuable information about the genetic diversity of the population and its evolutionary history.

SwarmGenomics compromises of a step-by-step process of whole-genome assembly, followed by modules which guide the user through different genomic analyses. At present, the modules are heterozygosity analysis, runs of homozygosity detection, PSMC analysis, unmapped reads analysis, repeat analysis, mitochondrial genome assembly, and nuclear mitochondrial DNA (NUMT) analysis. The modules include key population genetic and evolutionary analyses, from which the user can choose from. The pipeline has been created especially for researcher with limited bioinformatics experience and students, and can be used as a teaching material or as self-learning.

For more information about SwarmGenomics, please see our preprint (https://www.biorxiv.org/content/10.1101/2025.08.13.670070v1)



## Overview of the Pipeline
<img width="5113" height="2625" alt="swarmgenomics_pipeline" src="https://github.com/user-attachments/assets/f1b051a2-4cdf-4a0f-87c7-a97968415f1f" />


## Running SwarmGenomics with de.NBI, PuTTY, and FileZilla
One option for running SwarmGenomics is through the de.NBI Cloud (https://cloud.denbi.de/), which provides access to high-performance computing resources for bioinformatics projects. After setting up an account and requesting resources, you can connect to your virtual machine (VM) in different ways depending on your operating system:
- **On Windows**:
    - Use **PuTTY** to establish an SSH connection to the VM. This lets you log in to the remote machine, install packages, and run SwarmGenomics directly on the de.NBI infrastructure.
    - Use **FileZilla** to transfer files between your local computer and the VM, making it straightforward to upload scripts or download results.

- **On Linux or macOS**:
    - You can use the built-in Terminal for SSH access:
```
ssh -i /path/to/your/ssh/private/key ubuntu@123.456.78.90 -p 12345
```
## Get started

Start with the [Installation](0.%20Installations.md) instruction and [prepare your data](01.%20Getting%20started.md). Then follow the course description in sequential order.


## Input Requirements for SwarmGenomics Modules

The full SwarmGenomics pipeline generally requires **a reference genome** and **raw sequencing reads** (e.g., SRA files from short-read sequencing). However, many modules can be run independently using different types of input data.

The table below shows which types of input data can be used for each SwarmGenomics module.
- X = input is directly usable
- (✓) = input can be used, but only after preprocessing or extra steps (e.g. BAM → VCF before running heterozygosity).

| Input \ Module       | Genome features | Heterozygosity | Runs of Homozygosity | Genome visualisation | PSMC | Mitogenome assembly | NUMT identification | Repeat analysis | Unmapped reads |
| -------------------- | --------------- | -------------- | -------------------- | -------------------- | ---- | ------------------- | ------------------- | --------------- | -------------- |
| **BAM**              | X               | (✓)            | (✓)                  | X                    | (✓)  |                     |                     |                 | X              |
| **VCF**              |                 | X              | X                    | X                    | X    |                     |                     |                 |                |
| **FASTQ**            |                 |                |                      |                      |      | X                   |                     |                 |                |
| **Reference genome** |                 |                |                      | X                    |      |                     | X                   | X               |                |
| **Mitogenome**       |                 |                |                      |                      |      |                     | X                   |                 |                |
