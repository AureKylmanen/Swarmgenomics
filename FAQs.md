# Frequently asked questions
This section addresses common questions and issues you might encounter while using SwarmGenomics. It’s designed to help you get started quickly, troubleshoot common problems, and better understand how to use the different modules. If you don’t find your question here, please open an issue on the GitHub repository or check the documentation for the specific tool.

### What input files do I need?

Running through whole SwarmGenomics requires a reference genome and raw sequencing reads as SRA file or FASTQ. Most modules require a BAM, VCF, or FASTQ file, depending on the analysis step. Each module’s documentation specifies the required input.

### Can I use my own BAM/VCF/FASTQ files?

Yes. You can use files you’ve generated previously, as long as they are formatted correctly. Alternatively, you can follow the SwarmGenomics pipeline to generate them.

### I don’t have access to a powerful computer. Can I still run this?

Yes. You can run the pipeline on the de.NBI Cloud or other HPC clusters. Instructions are provided for connecting via terminal (Linux/Mac) or using PuTTY/FileZilla (Windows).

### Can I run only one part of the pipeline?

Yes. Each module is independent and can be run separately if you already have the required input files.

### My Bash script isn’t running, what should I do?

If your script doesn’t execute, it might not have the proper permissions or could have Windows-style line endings. You can fix this by:
```
# Make the script executable
chmod +x script_name.sh

# Convert Windows line endings to Unix
dos2unix script_name.sh
```
### How do I edit a file on the command line?

You can use a text editor in the terminal. Two common options are:

**nano (simple, beginner-friendly)**
```
nano filename.txt
```
- Use the arrow keys to navigate.
- Press Ctrl + O to save and Ctrl + X to exit.

**vim (more advanced)**
```
vim filename.txt
```
- Press i to enter insert mode and edit the file.
- Press Esc, then type :wq to save and quit.

### "No such file or directory" error

This usually means the software cannot find the file at the path you provided. Double-check that the file exists, that you’re in the correct working directory, and that the file name matches exactly (case-sensitive on Linux/Mac). If your file is in another folder, be sure to include the full path.

### PSMC isn’t working or gives me an error when I run it. What should I do?

One common issue is using a fragmented or low-coverage genome, which can lead to poor consensus sequences and unreliable demographic estimates. This typically leads to a large .psmcfa file (>200 MB). Another common pitfall is incorrect input formatting (e.g., problems when converting VCF to FASTQ), so double-check your preprocessing steps.

### SRA file doesn't download properly

Large SRA files can take a long time to download, and connections may occasionally drop or timeout. To prevent interruptions, you can use a small Bash script with nohup so the download continues even if you close your terminal:
```
#!/bin/bash
# download_sra.sh
prefetch SRRXXXXXXX
```
Then run it with:
```
# make executable:
chmod +x download_sra.sh

# run:
nohup bash download_sra.sh &
```
