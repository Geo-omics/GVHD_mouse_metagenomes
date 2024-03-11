import os
import re
import snakemake.io
from glob import glob

configfile: "config.yaml"
report: "code/report/workflow.rst"

# Substitute bash $USER environment variable with actual user id, otherwise some steps fail
param_work_dir = config["work_dir"] #get working directory from config file
userid = env_var = os.environ['USER'] #get bash $USER variable
work_dir = param_work_dir.replace("$USER", userid) #sub $USER for actual username

shell.prefix('printf "Job executed on: ${{HOSTNAME}}\n" && printf "SLURM job id: ${{SLURM_JOB_ID}}\n\n"; ')

# Get import sample names
#metaG_samples = glob_wildcards("import/metagenomes/{sample}/").sample
#metaT_samples = glob_wildcards("import/metatranscriptomes/{sample}/").sample
#metabolome_samples = glob_wildcards("import/metabolomes/{sample}/").sample
#amplicon_samples = glob_wildcards("import/amplicons/{sample}/").sample

# Get sample names
#metaG_samples = glob_wildcards("data/omics/metagenomes/{sample}/reads").sample
metaG_samples = open("data/sample_data/sample_names").read().splitlines()
#metaT_samples = glob_wildcards("data/omics/metatranscriptomes/{sample}/").sample
#metabolome_samples = glob_wildcards("data/omics/metabolomes/{sample}/").sample
#amplicon_samples = glob_wildcards("data/omics/amplicons/{sample}/").sample

# rule assemble:
#     input: 
#         expand("data/omics/metagenomes/{sample}/assembly/metaspades/contigs.fasta",sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.renamed.fa",sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/reads/dedup_reads_fwd.fastq.gz",sample = metaG_samples)


# rule prodigal:
#     input: expand("data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa", sample = metaG_samples)

# rule test:
#     input: expand("data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",sample = metaG_samples)
#     output: "test.out"

# # rule import:
# #     input:
# #     output:
# #     resources: cpus=1, mem_mb=8000, time_min=2880, mem_gb = 8
# #     shell:
# #         """

# #         """

# # rule gzip_fastx:

# # rule unzip_fastx:


# # rule zip_fastqs:
# #     input:  SCRATCH + "/01RAW_fqs/{sample}"
# #     output: temp(SCRATCH + "/02ZIPPED_fqs/{sample}")
# #     params: outdir = SCRATCH + "/02ZIPPED_fqs/"
# #     run:

# #         if wildcards.sample.endswith('.fastq'):
# #             shell("echo gzip {input}")
# #             shell("echo mv {input}.gz {params.outdir}")
# #         else:
# #             shell("mv {input} {params.outdir}")


# rule clumpify:
#     priority: 1
#     input: 
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz"
#     output: 
#         done = touch("data/omics/metagenomes/{sample}/reads/done.touch")
#     params:
#         clumped_fwd_reads = "data/omics/metagenomes/{sample}/reads/clumped_raw_fwd_reads.fastq.gz",
#         clumped_rev_reads = "data/omics/metagenomes/{sample}/reads/clumped_raw_rev_reads.fastq.gz",
#         temp_dir="/tmp/akiledal/gvhd/clumpify/{sample}"
#     conda: "config/conda_yaml/main.yaml"
#     benchmark:
#         "benchmarks/clumpify/{sample}.txt"
#     log: "logs/clumpify/{sample}_initial.log"
#     resources: cpus=16, mem_mb = lambda wildcards, attempt: attempt * 160000,
#         time_min=2880, 
#         #partition = "largemem"
#     shell:
#         """
#         #BBtools use more memory than given, reduce amount given by 20% to stay within job specs.
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)

#         TMPDIR={params.temp_dir}

#         mkdir -p {params.temp_dir}

#         clumpify.sh \
#             -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.fwd_reads} \
#             in2={input.rev_reads} \
#             out1={params.clumped_fwd_reads} \
#             out2={params.clumped_rev_reads} \
#             groups=24 \
#             zl=9 pigz \
#             usetmpdir=t \
#             tmpdir={params.temp_dir} \
#             t={resources.cpus} \
#             1>{log} 2>&1

#         rm {input.fwd_reads} {input.rev_reads}
#         mv {params.clumped_fwd_reads} {input.fwd_reads}
#         mv {params.clumped_rev_reads} {input.rev_reads}
#         """

# rule fastqc_raw:
#     input:
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz"
#     output:
#         #fwd_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/fwd_fastqc.html",
#         #rev_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/rev_fastqc.html"
#         touch("data/omics/metagenomes/{sample}/reads/fastqc_raw/.done")
#     conda:
#           "config/conda_yaml/fastqc.yaml"
#     shell:
#         """
#         mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_raw
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_raw -t {resources.cpus} {input.fwd_reads}
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_raw -t {resources.cpus} {input.rev_reads}
#         """

# rule run_fastqc_raw:
#     input: expand("data/omics/metagenomes/{sample}/reads/fastqc_raw/.done", sample = metaG_samples)

# rule fastp:
#     priority: 2
#     input: 
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
#         clumped = "data/omics/metagenomes/{sample}/reads/done.touch"
#     output: 
#         tmp_fwd = temp("data/omics/metagenomes/{sample}/reads/temp_fastp_fwd_reads.fastq.gz"),
#         tmp_rev = temp("data/omics/metagenomes/{sample}/reads/temp_fastp_rev_reads.fastq.gz"),
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/fastp_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/fastp_rev_reads.fastq.gz",
#         html_dedup = "data/omics/metagenomes/{sample}/reads/qc/fastp_dedup.html",
#         json_dedup = "data/omics/metagenomes/{sample}/reads/qc/fastp_dedup.json",
#         html = "data/omics/metagenomes/{sample}/reads/qc/fastp.html",
#         json = "data/omics/metagenomes/{sample}/reads/qc/fastp.json"
#     conda: "config/conda_yaml/fastp.yaml"
#     benchmark:
#         "benchmarks/fastp/{sample}.txt"
#     log: "logs/fastp/{sample}_fastp.log"
#     resources: cpus=36, mem_mb = lambda wildcards, attempt: attempt * 60000,
#         time_min=2880
#     shell:
#         """
        
#         # First deduplicate
#         fastp \
#             -i {input.fwd_reads} -I {input.rev_reads} \
#             -o {output.tmp_fwd} -O {output.tmp_rev} \
#             -h {output.html_dedup} -j {output.json_dedup} \
#             --thread {resources.cpus} \
#             -z 9 \
#             --dedup \
#             --dup_calc_accuracy 6
        
#         # Then trim and remove adapters
#         fastp \
#             -i {output.tmp_fwd} -I {output.tmp_rev} \
#             -o {output.fwd_reads} -O {output.rev_reads} \
#             -h {output.html} -j {output.json} \
#             --thread {resources.cpus} \
#             -z 9 \
#             --length_required 50 \
#             --n_base_limit 5 \
#             --low_complexity_filter --complexity_threshold 7 \
#             --detect_adapter_for_pe \
#             --correction \
#             --cut_front \
#             --cut_tail \
#             --cut_window_size=4 \
#             --cut_mean_quality 20 \
#             --overrepresentation_analysis
#         """

# rule fastqc_fastp:
#     input:
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/fastp_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/fastp_rev_reads.fastq.gz",
#     output:
#         #fwd_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/fwd_fastqc.html",
#         #rev_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/rev_fastqc.html"
#         touch("data/omics/metagenomes/{sample}/reads/fastqc_fastp/.done")
#     conda:
#           "config/conda_yaml/fastqc.yaml"
#     shell:
#         """
#         mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_fastp
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_fastp -t {resources.cpus} {input.fwd_reads}
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_fastp -t {resources.cpus} {input.rev_reads}
#         """


# rule multiqc:
#     input: 
#         "data/omics/metagenomes/{sample}/reads/fastqc_fastp/.done",
#         "data/omics/metagenomes/{sample}/reads/fastqc_decontam/.done",
#         "data/omics/metagenomes/{sample}/reads/fastqc_raw/.done",
#         "data/omics/metagenomes/{sample}/reads/fastqc_teal_decon/.done"
#     output: 
#         multiqc_dir = directory("data/omics/metagenomes/{sample}/reads/qc/multiqc")
#     conda: "config/conda_yaml/multiqc.yaml"
#     benchmark:
#         "benchmarks/multiqc/{sample}.txt"
#     log: "logs/multiqc/{sample}.log"
#     resources: cpus=1, mem_mb = 20000, time_min=2880
#     shell:
#         """
#         multiqc --interactive -d data/omics/metagenomes/{wildcards.sample}/reads/fastqc_* -o {output.multiqc_dir}
#         """


# rule run_fastqc_fastp:
#     input: expand("data/omics/metagenomes/{sample}/reads/fastqc_fastp/.done", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/reads/qc/multiqc",sample = metaG_samples)


# rule deduplicate:
#     input: 
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
#         clumped = "data/omics/metagenomes/{sample}/reads/done.touch"
#     output: 
#         dedup_interleaved = temp("data/omics/metagenomes/{sample}/reads/dedup_interleaved.fastq.gz"),
#         dedup_reads_fwd = "data/omics/metagenomes/{sample}/reads/dedup_reads_fwd.fastq.gz",
#         dedup_reads_rev = "data/omics/metagenomes/{sample}/reads/dedup_reads_rev.fastq.gz"
#     conda: "config/conda_yaml/main.yaml"
#     benchmark:
#         "benchmarks/dedup/{sample}.txt"
#     log: "logs/dedup/{sample}_dedup.log"
#     # resources: cpus=36, time_min=2880,
#     #     mem_mb = lambda wildcards, attempt: attempt * 170000,
#     #     #partition = "largemem"
#     resources: 
#         partition = "largemem",
#         cpus = 24, 
#         time_min = 7200,
#         mem_mb = 1000000
#     shell:
#         """
#         #BBtools use more memory than given, amount given by 20% to stay within job specs.
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)

#         dedupe.sh \
#             -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.fwd_reads} \
#             in2={input.rev_reads} \
#             out={output.dedup_interleaved} \
#             t={resources.cpus} \
#             1>{log} 2>&1

#         # Dedup only outputs interleaved files, this just converts back to paired
#         reformat.sh in={output.dedup_interleaved} \
#             out1={output.dedup_reads_fwd} \
#             out2={output.dedup_reads_rev}
#         """

# rule run_dedup:
#     input: 
#         expand("data/omics/metagenomes/{sample}/reads/dedup_reads_fwd.fastq.gz", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/reads/trimmed_fwd_reads.fastq.gz", sample = metaG_samples)


# rule get_contaminants:
#     output: 
#         human_genome = "data/reference/contaminants/human.fa.gz",
#         spike_ins = "data/reference/contaminants/spike-ins.fa"
#     resources: cpus = 1
#     shell:
#         """
#         wget -O {output.human_genome} http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz
#         wget -O {output.spike_ins} https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis/raw/master/spike-ins.fa
#         """

# rule bb_index_human:
#     input:
#         human_genome = rules.get_contaminants.output.human_genome
#     output:
#         index = directory("data/reference/contaminants/human/ref"),
#         bbmap_index_path = directory("data/reference/contaminants/human")
#     params:
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/bbmap_index.log"
#     benchmark:
#         "benchmarks/bb_index.txt"
#     resources: cpus = 8, mem_mb = 80000
#     shell:
#         """        
#         bbmap.sh \
#             ref={input.human_genome} \
#             path={output.bbmap_index_path} \
#             t={resources.cpus} \
#             1>{log} 2>&1
#         """

# rule get_mouse_genome:
#     output:
#         mouse_genome = "data/reference/contaminants/mouse.fa.gz"
#     shell:
#         """
#         wget -O {output.mouse_genome} http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
#         """

# rule bb_index_mouse:
#     input:
#         mouse_genome = rules.get_mouse_genome.output.mouse_genome
#     output:
#         index = directory("data/reference/contaminants/mouse/ref"),
#         bbmap_index_path = directory("data/reference/contaminants/mouse")
#     params:
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/bbmap_index.log"
#     benchmark:
#         "benchmarks/bb_index.txt"
#     resources: cpus = 8, mem_mb = 50000
#     shell:
#         """        
#         bbmap.sh \
#             ref={input.mouse_genome} \
#             path={output.bbmap_index_path} \
#             t={resources.cpus} \
#             1>{log} 2>&1
#         """

# rule combine_contams_and_index:
#     input:
#         mouse_genome = rules.get_mouse_genome.output.mouse_genome,
#         human_genome = rules.get_contaminants.output.human_genome
#     output:
#         contaminants = "data/reference/contaminants/merged.fa.gz",
#         index = directory("data/reference/contaminants/merged/ref"),
#         bbmap_index_path = directory("data/reference/contaminants/merged")
#     #params:
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/bbmap_index.log"
#     benchmark:
#         "benchmarks/bb_index.txt"
#     resources: cpus = 36, mem_mb = 100000
#     shell:
#         """       

#         cat {input.mouse_genome} > {output.contaminants}
#         cat {input.human_genome} >> {output.contaminants}

#         bbmap.sh \
#             ref={output.contaminants} \
#             path={output.bbmap_index_path} \
#             t={resources.cpus} \
#             1>{log} 2>&1
#         """

# rule trim_and_remove_adapters:
#     input:
#         dedup_reads_fwd = rules.deduplicate.output.dedup_reads_fwd,
#         dedup_reads_rev = rules.deduplicate.output.dedup_reads_rev,
#         adapters = "data/reference/contaminants/adapters.fa"
#     output:
#         trimmed_fwd = temp("data/omics/metagenomes/{sample}/reads/trimmed_fwd_reads.fastq.gz"),
#         trimmed_rev = temp("data/omics/metagenomes/{sample}/reads/trimmed_rev_reads.fastq.gz")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/trim_and_remove_adapters/{sample}.log"
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     benchmark:
#         "benchmarks/trim_and_remove_adapters/{sample}.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 175000, time_min = 2880
#     shell:
#         """
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

#         echo "\n\n***Trimming reads and removing adapters***\n\n" >> {log}
        
#         bbduk.sh -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.dedup_reads_fwd} \
#             in2={input.dedup_reads_rev} \
#             out1={output.trimmed_fwd} \
#             out2={output.trimmed_rev} \
#             t={resources.cpus} \
#             minlen=50 \
#             qtrim=rl \
#             trimq=15 \
#             ref={input.adapters} \
#             path={params.bbmap_index_path} \
#             ktrim=r k=23 mink=11 hdist=1 \
#             1>{log} 2>&1
#         """

# rule remove_spike_ins:
#     input:
#         spike_ins = rules.get_contaminants.output.spike_ins,
#         trimmed_fwd = rules.trim_and_remove_adapters.output.trimmed_fwd,
#         trimmed_rev = rules.trim_and_remove_adapters.output.trimmed_rev
#     output:
#         phix_rm_fwd = temp("data/omics/metagenomes/{sample}/reads/phix_fwd_reads.fastq.gz"),
#         phix_rm_rev = temp("data/omics/metagenomes/{sample}/reads/phix_rev_reads.fastq.gz")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/remove_spike_ins/{sample}.log"
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     benchmark:
#         "benchmarks/remove_spike_ins/{sample}.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 175000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

#         echo "\n\n***Removing spike-in***\n\n" >> {log}
        
#         bbduk.sh -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.trimmed_fwd} \
#             in2={input.trimmed_rev} \
#             outu1={output.phix_rm_fwd} \
#             outu2={output.phix_rm_rev} \
#             t={resources.cpus} k=31 hdist=1 \
#             ref={input.spike_ins} \
#             path={params.bbmap_index_path} \
#             1>>{log} 2>&1
#             """

# rule remove_contaminants:
#     input:
#         phix_rm_fwd = "data/omics/metagenomes/{sample}/reads/phix_fwd_reads.fastq.gz",
#         phix_rm_rev = "data/omics/metagenomes/{sample}/reads/phix_rev_reads.fastq.gz",
#         contaminants = rules.combine_contams_and_index.output.contaminants,
#         index_path = "data/reference/contaminants/merged"
#         #bbmap_index = "data/reference/contaminants/ref"
#     output:
#         decon_fwd = temp("data/omics/metagenomes/{sample}/reads/decon_fwd_reads.fastq.gz"),
#         decon_rev = temp("data/omics/metagenomes/{sample}/reads/decon_rev_reads.fastq.gz")
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/read_qc/{sample}.log"
#     benchmark:
#         "benchmarks/remove_contaminants/{sample}.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 175000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"
        
#         echo "\n\n***Removing contamination***\n\n" >> {log}
#         bbmap.sh \
#             in1={input.phix_rm_fwd} \
#             in2={input.phix_rm_rev} \
#             outu1={output.decon_fwd} \
#             outu2={output.decon_rev} \
#             t={resources.cpus} fast=t \
#             ref={input.contaminants} \
#             path={input.index_path} \
#             -Xmx${{bbmap_mem}}m -eoom \
#             1>>{log} 2>&1
#         """

# rule remove_spike_ins_fastp:
#     priority: 3
#     input:
#         spike_ins = rules.get_contaminants.output.spike_ins,
#         trimmed_fwd = rules.fastp.output.fwd_reads,
#         trimmed_rev = rules.fastp.output.rev_reads
#     output:
#         phix_rm_fwd = temp("data/omics/metagenomes/{sample}/reads/phix_fwd_reads_fastp.fastq.gz"),
#         phix_rm_rev = temp("data/omics/metagenomes/{sample}/reads/phix_rev_reads_fastp.fastq.gz")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/remove_spike_ins/{sample}_fastp.log"
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     benchmark:
#         "benchmarks/remove_spike_ins_fastp/{sample}.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 175000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

#         echo "\n\n***Removing spike-in***\n\n" >> {log}
        
#         bbduk.sh -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.trimmed_fwd} \
#             in2={input.trimmed_rev} \
#             outu1={output.phix_rm_fwd} \
#             outu2={output.phix_rm_rev} \
#             t={resources.cpus} k=31 hdist=1 \
#             ref={input.spike_ins} \
#             path={params.bbmap_index_path} \
#             1>>{log} 2>&1

#         echo "\n\n*** DONE ***\n\n" >> {log}
#         """

# rule remove_contaminants_fastp:
#     priority: 4
#     input:
#         phix_rm_fwd = rules.remove_spike_ins_fastp.output.phix_rm_fwd,
#         phix_rm_rev = rules.remove_spike_ins_fastp.output.phix_rm_rev,
#         contaminants = rules.combine_contams_and_index.output.contaminants,
#         index_path = "data/reference/contaminants/merged"
#         #bbmap_index = "data/reference/contaminants/ref"
#     output:
#         decon_fwd = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         decon_rev = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/read_qc/{sample}_fastp.log"
#     benchmark:
#         "benchmarks/remove_contaminants_fastp/{sample}.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 175000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"
        
#         echo "\n\n***Removing contamination***\n\n" >> {log}
#         bbmap.sh \
#             in1={input.phix_rm_fwd} \
#             in2={input.phix_rm_rev} \
#             outu1={output.decon_fwd} \
#             outu2={output.decon_rev} \
#             t={resources.cpus} \
#             fast=t \
#             zl=9 \
#             ref={input.contaminants} \
#             path={input.index_path} \
#             -Xmx${{bbmap_mem}}m -eoom \
#             1>>{log} 2>&1
#         """

# rule fastqc_decontam:
#     input:
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
#     output:
#         touch("data/omics/metagenomes/{sample}/reads/fastqc_decontam/.done")
#     conda:
#           "config/conda_yaml/fastqc.yaml"
#     shell:
#         """
#         mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_decontam
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_decontam -t {resources.cpus} {input.fwd_reads}
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_decontam -t {resources.cpus} {input.rev_reads}
#         """

# rule run_fastqc_decontam:
#     input: expand("data/omics/metagenomes/{sample}/reads/fastqc_decontam/.done", sample = metaG_samples)


# rule count_reads_fastp:
#     priority: 5
#     input:
#         raw_reads_fwd = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         raw_reads_rev = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
#         deduped_reads_fwd = "data/omics/metagenomes/{sample}/reads/temp_fastp_fwd_reads.fastq.gz",
#         deduped_reads_rev = "data/omics/metagenomes/{sample}/reads/temp_fastp_rev_reads.fastq.gz",
#         qual_filt_and_trimmed_fwd = rules.fastp.output.fwd_reads,
#         qual_filt_and_trimmed_rev = rules.fastp.output.rev_reads,
#         decon_reads_fwd = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         decon_reads_rev = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
#         bbnorm_fwd = "data/omics/metagenomes/{sample}/reads/bbnorm_fwd_reads.fastq.gz",
#         bbnorm_rev = "data/omics/metagenomes/{sample}/reads/bbnorm_rev_reads.fastq.gz"
#     output:
#         "data/omics/metagenomes/{sample}/reads/{sample}_read_count_fastp.tsv"
#     shell:
#         """
#         printf "read_state\tfwd_read_count\trev_read_count\n" > {output}
#         printf "raw_reads\t$(($(zcat {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "deduped_reads\t$(($(zcat {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "filt_and_trimmed_reads\t$(($(zcat {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(zcat {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "decon_reads\t$(($(zcat {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "bbnorm_reads\t$(($(zcat {input.bbnorm_fwd} | wc -l) / 4 ))\t$(($(zcat {input.bbnorm_rev} | wc -l) / 4 ))\n" >> {output}
#         """

# rule run_count_reads:
#     input: 
#         expand("data/omics/metagenomes/{sample}/reads/{sample}_read_count_fastp.tsv", sample=metaG_samples),
#         #expand("data/omics/metagenomes/{sample}/reads/{sample}_read_count.tsv", sample=metaG_samples),

# rule remove_poly_pairs:
#     input:
#         decon_fwd = ancient("data/omics/metagenomes/{sample}/reads/decon_fwd_reads.fastq.gz"),
#         decon_rev = ancient("data/omics/metagenomes/{sample}/reads/decon_rev_reads.fastq.gz")
#     output:
#         cleaned_fwd = "data/omics/metagenomes/{sample}/reads/cleaned_fwd_reads.fastq.gz",
#         cleaned_rev = "data/omics/metagenomes/{sample}/reads/cleaned_rev_reads.fastq.gz"
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/remove_poly_pairs/{sample}.log"
#     benchmark:
#         "benchmarks/remove_poly_pairs/{sample}.txt"
#     resources: cpus = 1, mem_mb = lambda wildcards, attempt: attempt * 175000, time_min = 2880
#     shell:
#         """
#         echo "\n\n***Running RemovePolyPairs.pl***\n\n" > {log}
        
#         perl code/RemovePolyPairs.pl {input.decon_fwd} {input.decon_rev} 50 {output.cleaned_fwd} {output.cleaned_rev} >> {log}
        
#         echo "\n\n***DONE***\n\n" >> {log}
#         """

# rule fastqc_Teal_decon:
#     input:
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/cleaned_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/cleaned_rev_reads.fastq.gz"
#     output:
#         touch("data/omics/metagenomes/{sample}/reads/fastqc_teal_decon/.done")
#     conda:
#           "config/conda_yaml/fastqc.yaml"
#     shell:
#         """
#         mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_teal_decon
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_teal_decon -t {resources.cpus} {input.fwd_reads}
#         fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_teal_decon -t {resources.cpus} {input.rev_reads}
#         """

# rule decon_reads:
#     input: 
#         expand("data/omics/metagenomes/{sample}/reads/cleaned_fwd_reads.fastq.gz", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/reads/fastp_fwd_reads.fastq.gz",sample = metaG_samples)



# rule count_reads_Teal:
#     input:
#         raw_reads_fwd = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         raw_reads_rev = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
#         deduped_reads_fwd = "data/omics/metagenomes/{sample}/reads/dedup_reads_fwd.fastq.gz",
#         deduped_reads_rev = "data/omics/metagenomes/{sample}/reads/dedup_reads_rev.fastq.gz",
#         qual_filt_and_trimmed_fwd = "data/omics/metagenomes/{sample}/reads/trimmed_fwd_reads.fastq.gz",
#         qual_filt_and_trimmed_rev = "data/omics/metagenomes/{sample}/reads/trimmed_rev_reads.fastq.gz",
#         decon_reads_fwd = "data/omics/metagenomes/{sample}/reads/cleaned_fwd_reads.fastq.gz",
#         decon_reads_rev = "data/omics/metagenomes/{sample}/reads/cleaned_rev_reads.fastq.gz"
#     output:
#         "data/omics/metagenomes/{sample}/reads/{sample}_read_count.tsv"
#     shell:
#         """
#         printf "read_state\tfwd_read_count\trev_read_count\n" > {output}
#         printf "raw_reads\t$(($(zcat {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "deduped_reads\t$(($(zcat {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "filt_and_trimmed_reads\t$(($(zcat {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(zcat {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output}
#         printf "decon_reads\t$(($(zcat {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
#         """

# rule bbnorm:
#     input:
#         fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
#         rev_reads = rules.remove_contaminants_fastp.output.decon_rev
#     output:
#         fwd_norm = temp("data/omics/metagenomes/{sample}/reads/bbnorm_fwd_reads.fastq.gz"),
#         rev_norm = temp("data/omics/metagenomes/{sample}/reads/bbnorm_rev_reads.fastq.gz")
#     params: "target=100 mindepth=2 bits=16 prefilter"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/bbnorm/{sample}.log"
#     benchmark:
#         "benchmarks/bbnorm/{sample}.txt"
#     resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 120000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"
        
#         bbnorm.sh \
#             -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.fwd_reads} \
#             in2={input.rev_reads} \
#             out1={output.fwd_norm} \
#             out2={output.rev_norm} \
#             t={resources.cpus} \
#             {params} \
#             1>>{log} 2>&1
#         """


# rule assemble_metaspades:
#     input:
#         fwd_reads = rules.bbnorm.output.fwd_norm,
#         rev_reads = rules.bbnorm.output.rev_norm
#     output:
#         assembly_dir = directory("data/omics/metagenomes/{sample}/assembly/metaspades"),
#         contigs = "data/omics/metagenomes/{sample}/assembly/metaspades/contigs.fasta"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/assembly/metaspades/{sample}.log"
#     benchmark: "benchmarks/metaspades/{sample}.txt"
#     resources: cpus = 24, time_min=20160, mem_mb = 150000
#         #mem_mb = lambda wildcards, attempt: attempt * 100000
#     shell:
#         """
#         metaspades.py -t {resources.cpus} -1 {input.fwd_reads} -2 {input.rev_reads} -o {output.assembly_dir} > {log}
#         """

# rule assemble_megahit:
#     input:
#         fwd_reads = rules.bbnorm.output.fwd_norm,
#         rev_reads = rules.bbnorm.output.rev_norm
#     output:
#         contigs = "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa",
#         #touch("data/omics/metagenomes/{sample}/assembly/megahit/.done")
#     params:
#         assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/assembly/megahit/{sample}.log"
#     benchmark: "benchmarks/megahit/{sample}.txt"
#     resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
#     shell:
#         """
#         rm -r {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
#         megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {params.assembly_dir} > {log}
#         """
        

# rule bbnorm_COassembly:
#     input:
#         fwd_reads = expand("data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz", sample = metaG_samples),
#         rev_reads = expand("data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz", sample = metaG_samples)
#     output:
#         fwd_cat = temp("data/omics/metagenomes/coassembly/fwd_reads.fastq.gz"),
#         rev_cat = temp("data/omics/metagenomes/coassembly/rev_reads.fastq.gz"),
#         fwd_norm = temp("data/omics/metagenomes/coassembly/bbnorm_fwd_reads.fastq.gz"),
#         rev_norm = temp("data/omics/metagenomes/coassembly/bbnorm_rev_reads.fastq.gz")
#     params: "target=100 mindepth=2 bits=16 prefilter"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/bbnorm/coassembly.log"
#     benchmark:
#         "benchmarks/bbnorm/coassembly.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 170000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"
        
#         cat data/omics/metagenomes/*/reads/decon_fwd_reads_fastp.fastq.gz > {output.fwd_cat}
#         cat data/omics/metagenomes/*/reads/decon_rev_reads_fastp.fastq.gz > {output.rev_cat}

#         bbnorm.sh \
#             -Xmx${{bbmap_mem}}m -eoom \
#             in1={output.fwd_cat} \
#             in2={output.rev_cat} \
#             out1={output.fwd_norm} \
#             out2={output.rev_norm} \
#             t={resources.cpus} \
#             {params} \
#             1>>{log} 2>&1
#         """

# ruleorder: COassemble_megahit > assemble_megahit

# rule COassemble_megahit:
#     input:
#         fwd_reads = rules.bbnorm_COassembly.output.fwd_norm,
#         rev_reads = rules.bbnorm_COassembly.output.rev_norm
#     output:
#         contigs = "data/omics/metagenomes/coassembly/assembly/megahit/final.contigs.fa",
#         #touch("data/omics/metagenomes/{sample}/assembly/megahit/.done")
#     params:
#         assembly_dir = "data/omics/metagenomes/coassembly/megahit"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/assembly/megahit/coassembly.log"
#     benchmark: "benchmarks/megahit/coassembly.txt"
#     resources: cpus = 36, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 170000
#     shell:
#         """
#         rm -r {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
#         megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {params.assembly_dir} > {log}
#         """

# rule rename_megahit_contigs:
#     input: 
#         script = "code/rename_megahit_contigs.R",
#         #assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit",
#         contigs = "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa",
#         #assembly_done = "data/omics/metagenomes/{sample}/assembly/megahit/.done"
#     output:
#         contigs = "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.renamed.fa",
#         #done = touch("data/omics/metagenomes/{sample}/assembly/megahit/.contigs_renamed")
#     container: "docker://eandersk/r_microbiome"
#     resources: cpus = 1, time_min=200, mem_mb = 50000
#     shell:
#         """
#         cd /home/kiledal/scratch_gdick1/GVHD
        
#         pwd
        
#         ./{input.script} {input.contigs}
#         """

# rule merge_assemblies:
#     input:
#         metaspades_contigs = rules.assemble_metaspades.output.contigs,
#         megahit_contigs = rules.rename_megahit_contigs.output.contigs,
#         merge_contigs_script = "code/Strain-Level_Metagenome_Analysis/Merge_Contigs.pl"
#     output:
#         concat_contigs = temp("data/omics/metagenomes/{sample}/assembly/contigs_concat.fa"),
#         dedup1 = temp("data/omics/metagenomes/{sample}/assembly/{sample}_dedup1.fa"),
#         dedup2 = temp("data/omics/metagenomes/{sample}/assembly/{sample}_dedup2.fa"),
#         dedup3 = temp("data/omics/metagenomes/{sample}/assembly/{sample}_dedup3.fa"),
#         dedup4 = temp("data/omics/metagenomes/{sample}/assembly/{sample}_dedup4.fa"),
#         dedup4_dot = temp("data/omics/metagenomes/{sample}/assembly/{sample}_graph4.dot"),
#         dedup5 = temp("data/omics/metagenomes/{sample}/assembly/{sample}_MERGED_CONTIGS.fa"),
#         dedup6 = temp("data/omics/metagenomes/{sample}/assembly/{sample}_MCDD.fa")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/assembly/merge/{sample}_merge_assemblies.log"
#     benchmark: "benchmarks/assembly/merge/{sample}.txt"
#     resources: cpus = 8, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
#     shell:
#         """
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
#         cat {input.metaspades_contigs} {input.megahit_contigs} > {output.concat_contigs}
#         dedupe.sh -da -Xmx${{bbmap_mem}}m -eoom in={output.concat_contigs} out={output.dedup1} tuc mid=99 minscaf=200 rnc=f ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f sort=length absorbcontainment=t mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 1>{log} 2>&1
#         dedupe.sh -da -Xmx${{bbmap_mem}}m -eoom in={output.dedup1} out={output.dedup2} tuc mid=99 minscaf=200 rnc=t ngn=t fo c pc=t fmj=t rc=t cc=t fcc=t mst=f sort=length absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 1>>{log} 2>&1
#         dedupe.sh -da -Xmx${{bbmap_mem}}m -eoom in={output.dedup2} out={output.dedup3} tuc mid=99 minscaf=200 rnc=t ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f ordered=t absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 1>>{log} 2>&1
#         dedupe.sh -da -Xmx${{bbmap_mem}}m -eoom in={output.dedup3} out={output.dedup4} tuc mid=99 minscaf=200 rnc=f ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f ordered=t absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t dot={output.dedup4_dot} t={resources.cpus} 1>>{log} 2>&1
#         perl {input.merge_contigs_script} data/omics/metagenomes/{wildcards.sample}/assembly/{wildcards.sample} 99 1>>{log} 2>&1
#         dedupe.sh -da -Xmx${{bbmap_mem}}m -eoom in={output.dedup5} out={output.dedup6} t={resources.cpus} tuc mid=99 minscaf=200 overwrite=f 1>>{log} 2>&1
#         """

# rule get_MEC:
#     output: 
#         mec_dir = directory("code/MEC"),
#         mec_script = "code/MEC/src/mec.py"
#     shell:
#         """
#         rm -rf {output.mec_dir}
#         cd code
#         git clone https://github.com/bioinfomaticsCSU/MEC.git
#         chmod +x ../{output.mec_script}
#         """

# rule correct_contigs:
#     input:
#         assembly = rules.merge_assemblies.output.dedup6,
#         fwd_reads = rules.remove_poly_pairs.output.cleaned_fwd,
#         rev_reads = rules.remove_poly_pairs.output.cleaned_rev,
#         mec_script = "code/MEC/src/mec.py"
#     output:
#         read_mapping = temp("data/omics/metagenomes/{sample}/assembly/{sample}_cleaned_reads_unsorted.sam"),
#         read_mapping_sorted = temp("data/omics/metagenomes/{sample}/assembly/{sample}_cleaned_reads.bam")
#     params:
#         assembly_dir = "data/omics/metagenomes/{sample}/assembly"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/assembly/correct_contigs/{sample}_correct_contigs.log"
#     benchmark: "benchmarks/assembly/correct_contigs/{sample}.txt"
#     resources: cpus = 8, mem_mb = lambda wildcards, attempt: attempt * 100000, time_min = 2880
#     shell:
#         """ 
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
#         mem_per_thread=$(( ({resources.mem_mb} - 2000) / {resources.cpus}))

#         #ALIGN READS TO CONTIGS
#         bbmap.sh -da -Xmx${{bbmap_mem}}m -eoom ref={input.assembly} path={params.assembly_dir} t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 idfilter=0.999 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} out={output.read_mapping} 1>{log} 2>&1
#         samtools view -bShu {output.read_mapping} | samtools sort -m ${{mem_per_thread}}M -@ {resources.cpus} -o {output.read_mapping_sorted} 1>>{log} 2>&1
#         samtools index {output.read_mapping_sorted} 1>>{log} 2>&1
#         """

# rule MEC:
#     input:
#         assembly = rules.merge_assemblies.output.dedup6,
#         read_mapping_sorted = rules.correct_contigs.output.read_mapping_sorted,
#         mec_script = rules.get_MEC.output.mec_script
#     output:        
#         corrected_assembly = "data/omics/metagenomes/{sample}/assembly/{sample}_merged_and_corrected_assembly.fa"
#     conda: "config/conda_yaml/mec.yaml"
#     log: "logs/assembly/MEC/{sample}.log"
#     benchmark: "benchmarks/assembly/MEC/{sample}.txt"
#     resources: cpus = 1, mem_mb = lambda wildcards, attempt: attempt * 40000, time_min = 2880
#     shell:
#         """
#         python {input.mec_script} \
#             -bam {input.read_mapping_sorted} \
#             -i {input.assembly} \
#             -o {output.corrected_assembly} 1>>{log} 2>&1
#         """

# rule run_assemblies:
#     input: expand("data/omics/metagenomes/{sample}/assembly/{sample}_merged_and_corrected_assembly.fa", sample=metaG_samples)

# rule quast:
#     input:
#         combined_contigs = rules.MEC.output.corrected_assembly,
#         metaspades_contigs = rules.assemble_metaspades.output.contigs,
#         megahit_contigs = rules.rename_megahit_contigs.output.contigs
#     output: directory("data/omics/metagenomes/{sample}/assembly/quast")
#     log: "logs/assembly/quast/{sample}.log"
#     benchmark:
#         "benchmarks/assembly/quast/{sample}.txt"
#     conda:
#         "config/conda_yaml/quast.yaml"
#     resources:
#         cpus = 1, mem_mb = 20000
#     shell:
#         """
#         quast.py {input.combined_contigs} {input.metaspades_contigs} {input.megahit_contigs} -o {output}
#         """

# rule run_quast:
#     input: expand("data/omics/metagenomes/{sample}/assembly/quast", sample = metaG_samples)


# rule predict_genes_and_calc_abundance:
#     input:
#         assembly = rules.rename_megahit_contigs.output.contigs,
#         fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
#         rev_reads = rules.remove_contaminants_fastp.output.decon_rev
#     output:
#         proteins = "data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa",
#         genes = "data/omics/metagenomes/{sample}/genes/{sample}_GENES.fna",
#         reads_vs_genes_rpkm = "data/omics/metagenomes/{sample}/genes/{sample}_READSvsGENES.rpkm",
#         reads_vs_contigs_rpkm = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.rpkm",
#         reads_vs_assembly_sam_gz = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam.gz"
#     params:
#         reads_vs_assembly_sam = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/predict_genes_and_calc_abundance/{sample}_predict_genes_and_calc_abundance.log"
#     benchmark: "benchmarks/predict_genes_and_calc_abundance/{sample}.txt"
#     resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 64000, time_min = 2880
#     shell:
#         """
#         prodigal -p meta -i {input.assembly} -a {output.proteins} -d {output.genes} 1>{log} 2>&1
#         bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {output.genes}) ref={output.genes} rpkm={output.reads_vs_genes_rpkm} 1>>{log} 2>&1
#         bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.assembly}) ref={input.assembly} rpkm={output.reads_vs_contigs_rpkm} 32bit=t outm={params.reads_vs_assembly_sam} 1>>{log} 2>&1
#         gzip {params.reads_vs_assembly_sam}
#         """

# rule download_uniref:
#     output: 
#         uniref100="data/reference/uniref/uniref100.fasta.gz"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/make_diamond_uniref_db/download_uniref.log"
#     resources: cpus = 1, mem_mb=1000
#     shell:
#         """
#         #cd data/reference/uniref
#         #wget -N https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
        
#         #mkdir data/reference/uniref
#         cp data/reference/uniref100.fasta.gz data/reference/uniref/
#         """

# rule make_diamond_uniref_db:
#     input: 
#         uniref100="data/reference/uniref/uniref100.fasta.gz"
#     output: 
#         uniref100_diamond="data/reference/uniref/uniref100.dmnd"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/make_diamond_uniref_db/make_diamond_uniref_db.log"
#     benchmark: "benchmarks/make_diamond_uniref_db/make_diamond_uniref_db.txt"
#     resources: cpus = 32, mem_mb=150000, time_min = 2880
#     shell:
#         """
#         diamond makedb --threads {resources.cpus} --in {input.uniref100} -d data/reference/uniref/uniref100
#         """

# rule align_to_uniref:
#     input:
#         diamond_db = rules.make_diamond_uniref_db.output.uniref100_diamond,
#         genes = rules.predict_genes_and_calc_abundance.output.genes
#     output:
#         gene_uniref_alignment = "data/omics/metagenomes/{sample}/{sample}_GENES.m8"
#     params:
#         "--top 0.5 --threads 10 --query-cover 50 --strand both -f 6 qseqid qlen sseqid slen qstart qend sstart send evalue pident mismatch qcovhsp scovhsp"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/align_to_uniref/{sample}_align_to_uniref.log"
#     benchmark: "benchmarks/align_to_uniref/{sample}.txt"
#     resources: cpus = 16, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 100000
#     shell:
#         """
#         diamond blastx \
#             -d {input.diamond_db} \
#             -q {input.genes} \
#             -o {output.gene_uniref_alignment} \
#             {params} 2>&1 | tee {log}
#         """

# rule make_uniref_alignments:
#     input: expand("data/omics/metagenomes/{sample}/{sample}_GENES.m8",sample = metaG_samples)

# # Required input for the annotation script: 
# # print "\t1. The gene or protein sequences:                      [sample]_GENES.fna\n";
# # print "\t2. The alignment file of the sequences to UMRAD:       [sample]_GENES.m8 \n";
# # print "\t3. The rpkm file for reads aligned to the genes:       [sample]_READSvsGENES.rpkm\n";
# # print "\t4. The rpkm file for reads aligned to the contigs:     [sample]_READSvsCONTIGS.rpkm\n";
# # print "\t5. The sample contigs sequences:                       [sample]_MCDD.fa\n\n";


# #             /TAXONOMY\_DB.*$year\.txt/){	$intax	=$refdir.$file;}
# # if($file =~ /UNIREF100\_INFO.*$year\.txt/){	$ininfo	=$refdir.$file;}
# # if($file =~ /Function\_Names.*\.txt/){	


# rule perl_dependencies:
#     output: touch("data/reference/perl_deps.touch")
#     shell: 
#         """
#         perl -MCPAN -Mlocal::lib -e 'CPAN::install(Number::Format)'
#         perl -MCPAN -Mlocal::lib -e 'CPAN::install(Statistics::Basic)'
#         perl -MCPAN -Mlocal::lib -e 'CPAN::install(Sort::Naturally)'
#         perl -MCPAN -Mlocal::lib -e 'CPAN::install(warnings)'

#         #export PERL5LIB=$CONDA_PREFIX/lib/perl5
#         #cpan CPAN::DistnameInfo
#         #cpan Number::Format
#         #cpan Statistics::Basic
#         #cpan Sort::Naturally
#         """
    
# rule annotate_contigs:
#     input:
#         genes = rules.predict_genes_and_calc_abundance.output.genes,
#         gene_alignment = rules.align_to_uniref.output.gene_uniref_alignment,
#         gene_rpkm = rules.predict_genes_and_calc_abundance.output.reads_vs_genes_rpkm,
#         contig_rpkm = rules.predict_genes_and_calc_abundance.output.reads_vs_contigs_rpkm,
#         contigs = rules.MEC.output.corrected_assembly,
#         script = "code/AnnotateContigs.pl",
#         UMRAD = "data/reference/UMRAD",
#         deps = "data/reference/perl_deps.touch"
#     output:
#         annotation_dir = directory("data/omics/metagenomes/{sample}/annotation")
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus=1, time_min = 2880, partition = "standard", mem_mb = lambda wildcards, attempt: attempt * 170000
#     ## for_debug  #resources: cpus=1, time_min = 30, partition = "debug", mem_mb = 8000 #mem_mb = lambda wildcards, attempt: attempt * 100000
#     log: "logs/annotate_contigs/{sample}.log"
#     benchmark: "benchmarks/annotate_contigs/{sample}.txt"
#     shell:
#         """
#         mkdir -p {output.annotation_dir}
#         ln {input.genes} {input.gene_alignment} {input.gene_rpkm} {input.contig_rpkm} {input.contigs} {output.annotation_dir}/
#         mv  {output.annotation_dir}/$(basename {input.contigs}) {output.annotation_dir}/{wildcards.sample}_MCDD.fa
#         proj_dir=$(pwd)

#         cp -f {input.script} {output.annotation_dir}/
#         ln {input.UMRAD}/* {output.annotation_dir}/

#         cd {output.annotation_dir}
#         perl AnnotateContigs.pl -s {wildcards.sample} #1>${{proj_dir}}/{log} 2>&1
        
#         #-d ${{proj_dir}}/{input.UMRAD}
#         #> ${{proj_dir}}/{log}
#         """

# rule teal_annotation:
#     input: 
#         #expand("results/metacodeR/{sample}_kraken.pdf", sample=metaG_samples),
#         expand("data/omics/metagenomes/{sample}/annotation", sample = metaG_samples),
#         #expand("data/omics/metagenomes/{sample}/bins", sample=metaG_samples)
    
# rule run_kraken:
#     input: expand("results/metacodeR/{sample}_kraken.pdf", sample=metaG_samples)

# rule initial_binning:
#     input:
#     output:
#     conda: "config/conda_yaml/main.yaml"
#     resources: 
#     shell:
#         """
#        while read i; do echo "doing $i"; 
#             #metabat
#             awk -F'\t' '{print $1"\t"$2"\t"$6"\t"$6"\t"$6"\t"}' ${i}_READSvsCONTIGS.rpkm > ${i}_MB_abund.txt
#             metabat2 --maxP=97 --minS=93 --maxEdges=300 -m 1500 -i ${i}_MERGED_CONTIGS_COR.fasta -a ${i}_MB_abund.txt -o ${i}_MB_bins
#             metabat2 --maxP=99 --minS=97 --maxEdges=300 -m 1500 -i ${i}_MERGED_CONTIGS_COR.fasta -a ${i}_MB_abund.txt -o ${i}_MB_bins
#             #maxbin
#             awk -F'\t' '{print $1"\t"$6}' ${i}_READSvsCONTIGS.rpkm | grep "CLUSTER" > ${i}_MX_abund.txt
#             perl run_MaxBin.pl -thread 20 -contig ${i}_MERGED_CONTIGS_COR.fasta -abund ${i}_MX_abund.txt -out ${i}_MX_bins
#         done < samp_list.txt;
#         """

# # Target rules for running kraken
# rule calc_kraken_uniq:
#     input: expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_out.txt", database = ["refseq","gtdb_r202"], sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_brackenMpa.txt", database = ["refseq","gtdb_r202"], sample = metaG_samples)

# rule calc_kraken_uniq_gtdb:
#     input: expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_brackenMpa.txt", database = ["gtdb_r202"], sample = metaG_samples)

# # Inspect Kraken2 databases, nescessary to produce nicely formatted taxonomy files for downstream analysis
# rule kraken_inspect:
#     input:
#         db = ancient("/geomicro/data2/kiledal/references/kraken_databases/{database}")
#     output: 
#         inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/{database}/inspect.txt"
#     conda: "code/main_env.yaml"
#     resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
#     shell:
#         """
#         kraken2-inspect --db {input.db} > {output.inspect_file}
#         """

# rule add_lineage_to_inspect_gtdb:
#     input:
#         db = ancient("/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202"),
#         inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect.txt"
#     output: 
#         inspect_w_lineage = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect_w_lineage.txt"
#     conda: "code/main_env.yaml"
#     resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
#     shell:
#         """
#         taxonkit lineage \
#             {input.inspect_file} \
#             --data-dir {input.db}/taxonomy \
#             -i 5 \
#             -o {output.inspect_w_lineage}
#         """

# rule add_lineage_to_inspect_refseq:
#     input:
#         db = ancient("/geomicro/data2/kiledal/references/kraken_databases/refseq"),
#         inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect.txt"
#     output: 
#         inspect_w_lineage_unformatted = temp("/geomicro/data2/kiledal/references/kraken_databases/refseq/unformatted_inspect_w_lineage.txt"),
#         inspect_w_lineage = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect_w_lineage.txt"
#     conda: "code/main_env.yaml"
#     resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
#     shell:
#         """
#         taxonkit lineage \
#             {input.inspect_file} \
#             --data-dir {input.db}/taxonomy \
#             -i 5 \
#             -o {output.inspect_w_lineage_unformatted}

#         taxonkit reformat \
#             {output.inspect_w_lineage_unformatted} \
#             -i 7 \
#             -P \
#             -o {output.inspect_w_lineage}
#         """


# rule kraken_database_tax_merge:
#     input:
#         script = "code/merge_kraken_tax.R",
#         gtdb_tax_info = rules.add_lineage_to_inspect_gtdb.output.inspect_w_lineage,
#         refseq_tax_info = rules.add_lineage_to_inspect_refseq.output.inspect_w_lineage
#     output:
#         combined_tax_info = "data/reference/kraken_tax_info_merged.tsv"
#     resources: cpus=1, mem_mb=4000, time_min=60
#     container: "docker://eandersk/r_microbiome"
#     shell:
#         """
#         ./{input.script} -g {input.gtdb_tax_info} -r {input.refseq_tax_info} -o {output.combined_tax_info}
#         """



# # Run kraken2 with KrakenUniq like functionality to screen for higher confidence results
# # First runs with a GTDB database to annotate bacterial and archaeal reads
# rule kraken2_gtdb_w_uniq:
#     input:
#         f_seq = rules.remove_poly_pairs.output.cleaned_fwd,
#         r_seq = rules.remove_poly_pairs.output.cleaned_rev,
#         # f_seq = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         # r_seq = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
#         db = "/scratch/gdick_root/gdick1/kiledal/kraken_databases/gtdb_r202",
#         kreport2mpa = "code/kreport2mpa.py"
#     output:
#         report = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_report.txt",
#         out = temp("data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_out.txt"),
#         bracken = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_bracken.txt",
#         bracken_report = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_brackenReport.txt",
#         bracken_mpa = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_brackenMpa.txt",
#         unclass_f = temp("data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_unclassified_1.fasta"),
#         unclass_r = temp("data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_unclassified_2.fasta"),
#         bracken_input = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_for_bracken.txt"
#     params:
#         uniq_minimizer_threshold = 150,
#         unclass_out = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_unclassified#.fasta"
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
#     shell:
#         """
#         kraken2 \
#             --threads {resources.cpus} \
#             --report-minimizer-data \
#             --report {output.report} \
#             --output {output.out} \
#             --db {input.db} \
#             --minimum-hit-groups 3 \
#             --unclassified-out {params.unclass_out} \
#             --paired {input.f_seq} {input.r_seq}

#         echo "Kraken complete, filtering..."

#         awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

#         echo "Filtering complete, running bracken..."

#         bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

#         ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages

#         echo "Bracken complete. Quitting."
#         """

# rule kraken2_refseq_w_uniq: ##Run kraken2
#     input:
#         f_seq = rules.kraken2_gtdb_w_uniq.output.unclass_f,
#         r_seq = rules.kraken2_gtdb_w_uniq.output.unclass_r,
#         db = "/scratch/gdick_root/gdick1/kiledal/kraken_databases/refseq",
#         kreport2mpa = "code/kreport2mpa.py"
#     output:
#         report = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_report.txt",
#         out = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_out.txt",
#         bracken = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_bracken.txt",
#         bracken_report = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_brackenReport.txt",
#         bracken_mpa = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_brackenMpa.txt",
#         bracken_input = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_for_bracken.txt"
#     params:
#         uniq_minimizer_threshold = 150
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
#     shell:
#         """
#         kraken2 \
#             --threads {resources.cpus} \
#             --report-minimizer-data \
#             --report {output.report} \
#             --output {output.out} \
#             --db {input.db} \
#             --minimum-hit-groups 3 \
#             --paired {input.f_seq} {input.r_seq}

#         awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

#         bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

#         ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages
#         """

# rule kraken2_gtdb_w_uniq_fastp:
#     input:
#         f_seq = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         r_seq = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
#         # f_seq = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         # r_seq = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
#         db = "/scratch/gdick_root/gdick1/kiledal/kraken_databases/gtdb_r202",
#         kreport2mpa = "code/kreport2mpa.py"
#     output:
#         report = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_report.txt",
#         out = temp("data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_out.txt"),
#         bracken = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_bracken.txt",
#         bracken_report = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_brackenReport.txt",
#         bracken_mpa = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_brackenMpa.txt",
#         unclass_f = temp("data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_unclassified_1.fasta"),
#         unclass_r = temp("data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_unclassified_2.fasta"),
#         bracken_input = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_for_bracken.txt"
#     params:
#         uniq_minimizer_threshold = 150,
#         unclass_out = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_unclassified#.fasta"
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
#     shell:
#         """
#         kraken2 \
#             --threads {resources.cpus} \
#             --report-minimizer-data \
#             --report {output.report} \
#             --output {output.out} \
#             --db {input.db} \
#             --minimum-hit-groups 3 \
#             --unclassified-out {params.unclass_out} \
#             --paired {input.f_seq} {input.r_seq}

#         echo "Kraken complete, filtering..."

#         awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

#         echo "Filtering complete, running bracken..."

#         bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

#         ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages

#         echo "Bracken complete. Quitting."
#         """


# # Any reads not annotated with the GTDB database are then annotated with a RefSeq database
# rule kraken2_refseq_w_uniq_fastp: ##Run kraken2
#     input:
#         f_seq = rules.kraken2_gtdb_w_uniq_fastp.output.unclass_f,
#         r_seq = rules.kraken2_gtdb_w_uniq_fastp.output.unclass_r,
#         db = "/scratch/gdick_root/gdick1/kiledal/kraken_databases/refseq",
#         kreport2mpa = "code/kreport2mpa.py"
#     output:
#         report = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_report.txt",
#         out = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_out.txt",
#         bracken = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_bracken.txt",
#         bracken_report = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_brackenReport.txt",
#         bracken_mpa = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_brackenMpa.txt",
#         bracken_input = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_for_bracken.txt"
#     params:
#         uniq_minimizer_threshold = 150
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
#     shell:
#         """
#         kraken2 \
#             --threads {resources.cpus} \
#             --report-minimizer-data \
#             --report {output.report} \
#             --output {output.out} \
#             --db {input.db} \
#             --minimum-hit-groups 3 \
#             --paired {input.f_seq} {input.r_seq}

#         awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

#         bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

#         ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages
#         """

# # Combine the kraken annotations and produce count table
# rule kraken_summarize:
#     input:
#         script = "code/merge_bracken.R",
#         kraken_results = expand("data/omics/metagenomes/{sample}/kraken_fastp/{database}_{sample}_bracken.txt", database = ["refseq","gtdb"], sample = metaG_samples),
#         combined_tax_info = rules.kraken_database_tax_merge.output.combined_tax_info
#     output:
#         counts = "data/sample_data/bracken_counts.tsv",
#         rel_abund = "data/sample_data/bracken_rel_abund.tsv"
#     resources: cpus=1, mem_mb=5000, time_min=60
#     container: "docker://eandersk/r_microbiome"
#     shell:
#         """
#         ./{input.script} --taxonomy={input.combined_tax_info} --counts-out={output.counts} --rel-out={output.rel_abund}
#         """


# # Target rule to make all the metacodeR plots
# rule plot_metacoders:
#     input: expand("results/metacodeR/{sample}_kraken.pdf", sample=metaG_samples)


# rule metacodeR:
#     input:
#         script = "code/plot_metacoder.R",
#         abund = "data/sample_data/bracken_rel_abund.tsv"
#         #metadata = "data/metadata.tsv"
#     params:
#         sample = "{sample}"
#     output: "results/metacodeR/{sample}_kraken.pdf"
#     resources: cpus=1, mem_mb=8000, time_min=60
#     container: "docker://eandersk/r_microbiome"
#     shell:
#         """
#         {input.script} --abund={input.abund} --sample={params.sample} --output={output}
#         """



# rule rename_for_bwa_mem:
#     input: 
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz"
#     output: 
#         fwd_reads = temp("data/omics/metagenomes/{sample}/reads/renamed_fwd_reads.fastq.gz"),
#         rev_reads = temp("data/omics/metagenomes/{sample}/reads/renamed_rev_reads.fastq.gz")
#     conda: "config/conda_yaml/main.yaml"
#     benchmark: "benchmarks/rename_reads/{sample}.txt"
#     log: "logs/rename_reads/{sample}.log"
#     resources: cpus=16, mem_mb = lambda wildcards, attempt: attempt * 50000,
#         time_min=400
#     shell:
#         """
#         #BBtools use more memory than given, reduce amount given by 20% to stay within job specs.
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)

#         rename.sh \
#             -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.fwd_reads} \
#             in2={input.rev_reads} \
#             prefix=x \
#             out1={output.fwd_reads} \
#             out2={output.rev_reads} \
#             zl=9 pigz \
#             t={resources.cpus} \
#             1>{log} 2>&1
#         """


# # rule make_bins:
# #     input:
# #         contigs = rules.MEC.output.corrected_assembly,
# #         fwd_reads = "data/omics/metagenomes/{sample}/reads/renamed_fwd_reads.fastq.gz",
# #         rev_reads = "data/omics/metagenomes/{sample}/reads/renamed_rev_reads.fastq.gz"
# #     output:
# #         #directory("data/bins/{short_sample}/metaspades/concoct_bins"),
# #         #directory("data/bins/{short_sample}/metaspades/maxbin2_bins"),
# #         #directory("data/bins/{short_sample}/metaspades/metabat2_bins"),
# #         directory("data/omics/metagenomes/{sample}/bins/work_files"),
# #         "data/omics/metagenomes/{sample}/bins/work_files/assembly.fa",
# #         fwd_reads = temp("data/omics/metagenomes/{sample}/reads/reads_1.fastq"),
# #         rev_reads = temp("data/omics/metagenomes/{sample}/reads/reads_2.fastq"),
# #         out_dir = directory("data/omics/metagenomes/{sample}/bins")
# #     params:
# #         #out_dir = "data/omics/metagenomes/{sample}/bins"
# #     conda: "code/metawrap.yaml"
# #     log: "logs/binning/{sample}/binning.log"
# #     benchmark: "benchmarks/binning/{sample}/binning_benchmark.log"
# #     resources: cpus=36, mem_mb=50000, time_min=4320, mem_gb = 50
# #     shell:
# #         """
# #         WORK_DIR=$PWD

# #         #data/omics/metagenomes/{wildcards.sample}/bins
# #         #mkdir -p data/omics/metagenomes/{wildcards.sample}/bins
# #         #cp data/qc_sequence_files/{wildcards.sample}[a-z]_R[12].fastq.gz data/bins/{wildcards.sample}/metaspades/
# #         gunzip -c {input.fwd_reads} > {output.fwd_reads}
# #         gunzip -c {input.rev_reads} > {output.rev_reads}
 
# #         #cd data/bins/{wildcards.sample}/metaspades
# #         #rename _R1.fastq _1.fastq *_R1.fastq
# #         #rename _R2.fastq _2.fastq *_R2.fastq

# #         #cd $WORK_DIR

# #         metawrap binning \
# #             -a {input.contigs} \
# #             -o {output.out_dir} \
# #             -t {resources.cpus} \
# #             -m {resources.mem_gb} \
# #             --metabat2 \
# #             --maxbin2 \
# #             --concoct \
# #             --universal \
# #             {output.fwd_reads} \
# #             {output.rev_reads} 1>{log} 2>&1

# #         #rm -f data/bins/{wildcards.sample}/metaspades/*.fastq
# #         """

# rule sample_annotation:
#     input: 
#         expand("data/omics/metagenomes/{sample}/annotation", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/metabat2_bins", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/gtdbtk", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/checkm.txt", sample = metaG_samples)


# # rule initial_binning:
# #     input:
# #     output:
# #     conda: "config/conda_yaml/main.yaml"
# #     resources: 
# #     shell:
# #         """
# #        while read i; do echo "doing $i"; 
# #             #metabat
# #             awk -F'\t' '{print $1"\t"$2"\t"$6"\t"$6"\t"$6"\t"}' ${i}_READSvsCONTIGS.rpkm > ${i}_MB_abund.txt
# #             metabat2 --maxP=97 --minS=93 --maxEdges=300 -m 1500 -i ${i}_MERGED_CONTIGS_COR.fasta -a ${i}_MB_abund.txt -o ${i}_MB_bins
# #             metabat2 --maxP=99 --minS=97 --maxEdges=300 -m 1500 -i ${i}_MERGED_CONTIGS_COR.fasta -a ${i}_MB_abund.txt -o ${i}_MB_bins
# #             #maxbin
# #             awk -F'\t' '{print $1"\t"$6}' ${i}_READSvsCONTIGS.rpkm | grep "CLUSTER" > ${i}_MX_abund.txt
# #             perl run_MaxBin.pl -thread 20 -contig ${i}_MERGED_CONTIGS_COR.fasta -abund ${i}_MX_abund.txt -out ${i}_MX_bins
# #         done < samp_list.txt;
# #         """



# rule make_bins:
#     input:
#         contigs = rules.MEC.output.corrected_assembly,
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/renamed_fwd_reads.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/renamed_rev_reads.fastq.gz"
#     output:
#         directory("data/omics/metagenomes/{sample}/bins/concoct_bins"),
#         directory("data/omics/metagenomes/{sample}/bins/maxbin2_bins"),
#         directory("data/omics/metagenomes/{sample}/bins/metabat2_bins"),
#         directory("data/omics/metagenomes/{sample}/bins/work_files"),
#         "data/omics/metagenomes/{sample}/bins/work_files/assembly.fa",
#         bam = temp("data/omics/metagenomes/{sample}/bins/work_files/reads.bam"),
#         fwd_reads = temp("data/omics/metagenomes/{sample}/reads/reads_1.fastq"),
#         rev_reads = temp("data/omics/metagenomes/{sample}/reads/reads_2.fastq"),
#         #out_dir = directory("data/omics/metagenomes/{sample}/bins")
#     params:
#         out_dir = "data/omics/metagenomes/{sample}/bins"
#     conda: "config/conda_yaml/metawrap.yaml"
#     resources: cpus=16, mem_mb=50000, time_min=2880, mem_gb = 50
#     shell:
#         """
#         WORK_DIR=$PWD

#         #data/omics/metagenomes/{wildcards.sample}/bins
#         #mkdir -p data/omics/metagenomes/{wildcards.sample}/bins
#         #cp data/qc_sequence_files/{wildcards.sample}[a-z]_R[12].fastq.gz data/bins/{wildcards.sample}/metaspades/
#         gunzip -c {input.fwd_reads} > {output.fwd_reads}
#         gunzip -c {input.rev_reads} > {output.rev_reads}
 
#         #cd data/bins/{wildcards.sample}/metaspades
#         #rename _R1.fastq _1.fastq *_R1.fastq
#         #rename _R2.fastq _2.fastq *_R2.fastq

#         #cd $WORK_DIR

#         metawrap binning \
#             -a {input.contigs} \
#             -o {params.out_dir} \
#             -t {resources.cpus} \
#             -m {resources.mem_gb} \
#             --metabat2 \
#             --maxbin2 \
#             --concoct \
#             --universal \
#             {output.fwd_reads} \
#             {output.rev_reads}

#         #rm -f data/bins/{wildcards.sample}/metaspades/*.fastq
#         """


# rule link_reads_w_sample_names:
#     input: "data/omics/metagenomes/{sample}/reads/decon_{dir}_reads_fastp.fastq.gz"
#     output: "data/omics/metagenomes/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz"
#     resources: cpus=1, mem_mb = 500
#     shell:
#         """
#         ln {input} {output}
#         """


# rule map_to_contigs:
#     input:
#         expand("data/omics/metagenomes/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz", sample = metaG_samples, dir = ["fwd", "rev"]),
#         contigs = rules.rename_megahit_contigs.output.contigs
#     output: 
#         bam_dir = directory(temp("data/omics/metagenomes/{sample}/bins/bam"))
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=16, mem_mb=50000, time_min=2880, mem_gb = 50
#     shell:
#         """
#         coverm make -c data/omics/metagenomes/*/reads/fastp_decon/*.fastq.gz \
#             -r {input.contigs} \
#             --discard-unmapped \
#             -t {resources.cpus} \
#             -o {output.bam_dir} 
#         """

# rule run_contig_map:
#     input: expand("data/omics/metagenomes/{sample}/bins/bam", sample = metaG_samples)

# rule contig_coverage:
#     input:
#         bam_dir = "data/omics/metagenomes/{sample}/bins/bam"
#     output: 
#         coverage = "data/omics/metagenomes/{sample}/bins/contig_coverage.tsv",
#         coverage_metabat = "data/omics/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=8, mem_mb=100000, time_min=2880
#     shell:
#         """
#         coverm contig \
#             -b {input.bam_dir}/*.bam \
#             -t {resources.cpus} \
#             --output-file {output.coverage}

#         coverm contig \
#             -b {input.bam_dir}/*.bam \
#             -t {resources.cpus} \
#             --methods metabat \
#             --output-file {output.coverage_metabat}
#         """


# rule index_contig_coverage:
#     input:
#         bam_dir = "data/omics/metagenomes/{sample}/bins/bam"
#     output: 
#         index_done = touch("data/omics/metagenomes/{sample}/bins/.bam_indexed")
#     params:
#         bam = "data/omics/metagenomes/{sample}/bins/bam/final.contigs.renamed.fa.decon_fwd_reads_fastp.fastq.gz.bam"
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=4, mem_mb=20000, time_min=2880
#     shell:
#         """
#         parallel -j {resources.cpus} samtools index -@ 1 ::: {input.bam_dir}/*.bam
#         """

# rule concoct:
#     input:
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         bam_index = rules.index_contig_coverage.output.index_done,
#         bam_dir = "data/omics/metagenomes/{sample}/bins/bam"
#     output:
#         cut_contigs = "data/omics/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa",
#         cut_contigs_bed = "data/omics/metagenomes/{sample}/bins/CONCOCT/contigs_10K.bed",
#         cut_coverage = "data/omics/metagenomes/{sample}/bins/CONCOCT/cut_coverage_table.tsv"
#     params:
#         outdir = "data/omics/metagenomes/{sample}/bins/CONCOCT/output",
#         bam = "data/omics/metagenomes/{sample}/bins/bam/*.bam"
#     benchmark: "benchmarks/concoct/{sample}.txt"
#     conda: "config/conda_yaml/concoct.yaml"
#     resources: cpus=16, mem_mb=100000, time_min=10080, mem_gb = 50
#     shell:
#         """
#         cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output.cut_contigs_bed} > {output.cut_contigs}
        
#         concoct_coverage_table.py {output.cut_contigs_bed} {params.bam} > {output.cut_coverage}

#         concoct --threads {resources.cpus} --composition_file {output.cut_contigs} --coverage_file {output.cut_coverage} -b {params.outdir}/
        
#         merge_cutup_clustering.py {params.outdir}/clustering_gt1000.csv > {params.outdir}/clustering_merged.csv

#         mkdir -p {params.outdir}/fasta_bins
#         extract_fasta_bins.py {input.contigs} {params.outdir}/clustering_merged.csv --output_path {params.outdir}/fasta_bins
#         """


# rule metabat2:
#     input:
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         bam_index = rules.index_contig_coverage.output.index_done,
#         bam_dir = "data/omics/metagenomes/{sample}/bins/bam",
#         coverm_depth = "data/omics/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
#     output:
#         #depth = "data/omics/metagenomes/{sample}/bins/jgi_depth_summary.txt",
#         done = touch("data/omics/metagenomes/{sample}/bins/METABAT2/.done")
#     params:
#         bin_name = directory("data/omics/metagenomes/{sample}/bins/METABAT2/metabat2")
#     benchmark: "benchmarks/metabat2/{sample}.txt"
#     singularity: "docker://metabat/metabat"
#     resources: cpus=16, mem_mb=20000, time_min=2880, mem_gb = 50
#     shell:
#         """
#         pwd 

#         cd /home/kiledal/scratch_gdick1/GVHD
        
#         pwd

#         metabat2 -i {input.contigs} -a {input.coverm_depth} -o {params.bin_name} -m 2000 -t {resources.cpus} --unbinned
#         """

#         #jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam_dir}/*.bam

# rule maxbin2_coverage:
#     input:
#         script = "code/create_maxbin_coverage.R",
#         coverm_depth = "data/omics/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
#     output:
#         depths_file = "data/omics/metagenomes/{sample}/bins/maxbin/depths.txt"
#     singularity: "docker://eandersk/r_microbiome"
#     resources: cpus=1, mem_mb=50000, time_min=1000
#     shell:
#         """
#         cd /home/kiledal/scratch_gdick1/GVHD
        
#         pwd

#         ./{input.script} {input.coverm_depth}
#         """

# rule maxbin2:
#     input:
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         depth = rules.maxbin2_coverage.output.depths_file
#     output:
#         done = touch("data/omics/metagenomes/{sample}/bins/maxbin/.done")
#     params:
#         bin_dir = "data/omics/metagenomes/{sample}/bins/maxbin/maxbin"
#     benchmark: "benchmarks/maxbin/{sample}.txt"
#     conda: "config/conda_yaml/maxbin.yaml"
#     resources: cpus=16, mem_mb=20000, time_min=10080, mem_gb = 50
#     shell:
#         """
#         pwd 

#         cd /home/kiledal/scratch_gdick1/GVHD
        
#         pwd
        
#         run_MaxBin.pl -contig {input.contigs} \
#             -markerset 107 \
#             -thread {resources.cpus} \
#             -min_contig_length 2000\
# 	        -out {params.bin_dir} \
# 	        -abund_list {input.depth}
#         """

# #### For testing on only one sample ####
# #metaG_samples = "c39841b318d0487eda9e0134e5c06381"
# #metaG_samples = "coassembly"
# ########################################


# rule calc_contig_coverage:
#     input: 
#         expand("data/omics/metagenomes/{sample}/bins/contig_coverage.tsv", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.renamed.fa", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/semibin", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/METABAT2/.done/", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/maxbin/.done", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/VAMB", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/metadecoder/.done", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm.txt", sample =metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/das_tool/.done", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/.drep_done", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/.done_GTDB", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/.done_gunc", sample = metaG_samples)

# rule run_checkm_new_per_sample:
#     input: expand("data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm.txt", sample =metaG_samples)


# rule gunc_and_drep:
#     input: 
#         expand("data/omics/metagenomes/{sample}/bins/.done_gunc", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/bins/.drep_done", sample = metaG_samples)


# rule bin_coassembly:
#     input: 
#         expand("data/omics/metagenomes/{sample}/bins/contig_coverage.tsv", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.renamed.fa", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/semibin", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/METABAT2/.done/", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/maxbin/.done", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/VAMB", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/metadecoder/.done", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm.txt", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/das_tool/.done", sample = "coassembly"),
#         expand("data/omics/metagenomes/{sample}/bins/.drep_done", sample = "coassembly")



# rule semibin_ref_download:
#     output: directory("data/reference/semibin/gtdb")
#     params:
#     conda: "config/conda_yaml/semibin.yaml"
#     resources: cpus =1, mem_mb = 2000, time_min=2880
#     shell:
#         """
#         SemiBin download_GTDB --reference-db {output}
#         """


# rule semibin:
#     input:
#         #"data/reference/semibin/.gtdb_downloaded",
#         ref_db = "data/reference/semibin/gtdb",
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         bam_dir = "data/omics/metagenomes/{sample}/bins/bam"
#     output:
#         out_dir = directory("data/omics/metagenomes/{sample}/bins/semibin"),
#         done = touch("data/omics/metagenomes/{sample}/bins/semibin/.done")
#     params:
#     conda: "config/conda_yaml/semibin.yaml"
#     benchmark: "benchmarks/semibin/{sample}.txt"
#     log: "logs/semibin/{sample}.log"
#     resources: cpus=16, mem_mb=170000, time_min=2880, mem_gb = 50
#     shell:
#         """
#         WORK_DIR=$PWD

#         SemiBin \
#             single_easy_bin \
#             --threads {resources.cpus} \
#             --reference-db {input.ref_db} \
#             -i {input.contigs} \
#             -b {input.bam_dir}/*.bam \
#             -o {output.out_dir} | tee {log}

#             #  --environment mouse_gut ## can only be used for single sample binning
#         """

# rule run_semibin:
#     input: expand("data/omics/metagenomes/{sample}/bins/semibin/.done", sample = metaG_samples)


# rule VAMB:
#     input:
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         coverm_depth = "data/omics/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
#     output:
#         outdir = directory("data/omics/metagenomes/{sample}/bins/VAMB")
#     #conda: "config/conda_yaml/VAMB.yaml"
#     conda: "VAMB"
#     benchmark: "benchmarks/VAMB/{sample}.txt"
#     log: "logs/VAMB/{sample}.log"
#     resources: cpus=1, mem_mb=40000, time_min=1440, partition = "gpu", gpu = 1
#     shell:
#         """
#         vamb -o _ --outdir {output.outdir} --fasta {input.contigs} --jgi {input.coverm_depth} --minfasta 200000 --cuda
#         """

# rule format_coverage_for_metadecoder:
#     input:
#         script = "code/make_metadecoder_coverage.R",
#         coverage = "data/omics/metagenomes/{sample}/bins/contig_coverage.tsv",
#         contigs = rules.rename_megahit_contigs.output.contigs
#     output: "data/omics/metagenomes/{sample}/bins/metadecoder/coverage.tsv"
#     singularity: "docker://eandersk/r_microbiome"
#     resources: cpus=1, mem_mb = 50000, time_min=360
#     shell:
#         """
#         pwd && cd /home/kiledal/scratch_gdick1/GVHD && pwd

#         ./{input.script} --coverage={input.coverage} --contigs={input.contigs} --out={output}
#         """
        
# rule metadecoder:
#     input:
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         coverage = "data/omics/metagenomes/{sample}/bins/metadecoder/coverage.tsv"
#     output:
#         touch("data/omics/metagenomes/{sample}/bins/metadecoder/.done"),
#         seed = "data/omics/metagenomes/{sample}/bins/metadecoder/seed.txt",
#         bins_dir = directory("data/omics/metagenomes/{sample}/bins/metadecoder/bins")
#     params:
#         out_prefix = "data/omics/metagenomes/{sample}/bins/metadecoder/bins/{sample}"
#     #conda: "config/conda_yaml/VAMB.yaml"
#     conda: "metadecoder"
#     shadow: "minimal"
#     benchmark: "benchmarks/metadecoder/{sample}.txt"
#     log: "logs/metadecoder/{sample}.log"
#     resources: cpus=1, mem_mb=150000, time_min=10080, partition = "gpu", gpu = 1
#     shell:
#         """
#         mkdir -p {output.bins_dir}
        
#         metadecoder seed --threads {resources.cpus} -f {input.contigs} -o {output.seed}

#         metadecoder cluster -f {input.contigs} -c {input.coverage} -s {output.seed}  -o {params.out_prefix} | tee {log}
#         """


# rule standardize_bins:
#     input:
#         "data/omics/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa",
#         "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.renamed.fa",
#         "data/omics/metagenomes/{sample}/bins/semibin",
#         "data/omics/metagenomes/{sample}/bins/METABAT2/.done/",
#         "data/omics/metagenomes/{sample}/bins/maxbin/.done", 
#         "data/omics/metagenomes/{sample}/bins/VAMB",
#         "data/omics/metagenomes/{sample}/bins/metadecoder/.done",
#         script = "code/standardize_bins.R"
#     output: 
#         contig_bin_mapping = "data/omics/metagenomes/{sample}/bins/contig_bins.rds",
#         bins_linked = "data/omics/metagenomes/{sample}/bins/all_raw_bins/.bins_linked"
#     params: 
#         sample = "{sample}"
#     singularity: "docker://eandersk/r_microbiome"
#     resources: cpus=1, mem_mb = 50000, time_min=360
#     shell:
#         """
#         pwd && cd /home/kiledal/scratch_gdick1/GVHD && pwd

#         ./{input.script} --sample={params.sample}
#         """

# rule checkm_new_per_sample:
#     input: "data/omics/metagenomes/{sample}/bins/all_raw_bins/.bins_linked"
#     output:
#         results = "data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm.txt"
#     params:
#         in_dir = "data/omics/metagenomes/{sample}/bins/all_raw_bins",
#         out_dir = "data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm"
#     conda: "config/conda_yaml/checkm.yaml"
#     resources: cpus=16, mem_mb=80000, time_min=2880
#     shell:
#         """
#         checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
#         """


# rule make_das_and_drep_inputs:
#     input:
#         "data/omics/metagenomes/{sample}/bins/all_raw_bins/.bins_linked",
#         "data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm.txt",
#         contig_bin_mapping = "data/omics/metagenomes/{sample}/bins/contig_bins.rds",
#         script = "code/make_das_and_drep_inputs.R"
#     output: 
#         metabat2_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/metabat2_contigs.tsv",
#         maxbin_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/maxbin_contigs.tsv",
#         concoct_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/concoct_contigs.tsv",
#         metadecoder_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/metadecoder_contigs.tsv",
#         semibin_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/semibin_contigs.tsv",
#         VAMB_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/VAMB_contigs.tsv",
#         drep_bin_info = "data/omics/metagenomes/{sample}/bins/bins_for_drep/genome_info.csv",
#         drep_bins_linked = "data/omics/metagenomes/{sample}/bins/bins_for_drep/.bins_linked"
#     params: 
#         sample = "{sample}"
#     singularity: "docker://eandersk/r_microbiome"
#     resources: cpus=1, mem_mb = 50000, time_min=360
#     shell:
#         """
#         pwd && cd /home/kiledal/scratch_gdick1/GVHD && pwd

#         ./{input.script} --sample={params.sample}
#         """


# rule checkm_new:
#     #input: "data/omics/metagenomes/metagenome_bins/raw_combined_bins"
#     output:
#         #dir = temp(directory("data/omics/metagenomes/metagenome_bins/raw_combined_bins")),
#         results = "data/omics/metagenomes/metagenome_bins/raw_combined_bins/checkm.txt"
#     params:
#         in_dir = "data/omics/metagenomes/metagenome_bins/raw_combined_bins",
#         out_dir = "data/omics/metagenomes/metagenome_bins/raw_combined_bins/checkm"
#     conda: "config/conda_yaml/checkm.yaml"
#     resources: cpus=24, mem_mb=120000, time_min=2880
#     shell:
#         """
#         checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
#         """


# rule dastool_new:
#     input:
#         contigs = rules.rename_megahit_contigs.output.contigs,
#         checkm_res = "data/omics/metagenomes/{sample}/bins/all_raw_bins/checkm.txt",
#         metabat2_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/metabat2_contigs.tsv",
#         maxbin_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/maxbin_contigs.tsv",
#         concoct_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/concoct_contigs.tsv",
#         metadecoder_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/metadecoder_contigs.tsv",
#         semibin_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/semibin_contigs.tsv",
#         VAMB_contigs = "data/omics/metagenomes/{sample}/bins/das_tool/VAMB_contigs.tsv"
#     params:
#         bin_folder = rules.make_bins.params.out_dir,
#         das_prefix = "data/omics/metagenomes/{sample}/bins/das_tool/output/{sample}"
#     output: 
#         #summary = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_summary.txt",
#         das_done = touch("data/omics/metagenomes/{sample}/bins/das_tool/.done"),
#         #das_bins_folder = directory("data/omics/metagenomes/{sample}/bins/das_tool/output/_DASTool_bins")
#     conda: "config/conda_yaml/das_tool.yaml"
#     benchmark: "benchmarks/dastool/{sample}.txt"
#     log: "logs/dastool/{sample}.log"
#     resources: cpus=8, mem_mb=50000, time_min=2880, mem_gb = 50
#     shell:
#         """
#         mkdir -p $(dirname {params.das_prefix})
        
#         DAS_Tool \
#             -i {input.metabat2_contigs},{input.maxbin_contigs},{input.concoct_contigs},{input.metadecoder_contigs},{input.semibin_contigs},{input.VAMB_contigs} \
#             -c {input.contigs} \
#             -o {params.das_prefix} \
#             -l metabat2,maxbin,concoct,metadecoder,semibin,VAMB \
#             --threads {resources.cpus} \
#             --write_bins \
#              | tee {log}
#         """


# rule drep_new:
#     input: 
#         "data/omics/metagenomes/{sample}/bins/bins_for_drep/.bins_linked"
#     output:
#         touch("data/omics/metagenomes/{sample}/bins/.drep_done")
#     params:
#         bins_linked = "data/omics/metagenomes/{sample}/bins/bins_for_drep/.bins_linked",
#         input_bins = "data/omics/metagenomes/{sample}/bins/bins_for_drep/*.fa",
#         genome_info = "data/omics/metagenomes/{sample}/bins/bins_for_drep/genome_info.csv", # Will need to be moved to inputs 
#         main_dir = directory("data/omics/metagenomes/{sample}/bins/drep/"),
#         MAGs= directory("data/omics/metagenomes/metagenome_bins/derep/dereplicated_genomes")
#     conda: "config/conda_yaml/drep.yaml"
#     resources: cpus=8, mem_mb=150000, time_min=2880
#     shell:
#         """
#         rm -rf {params.main_dir} # Clear any old drep output
        
#         dRep dereplicate \
#             {params.main_dir} \
#             -p {resources.cpus} \
#             --contamination 50 \
#             --completeness 30 \
#             --length 10000 \
#             --genomeInfo {params.genome_info} \
#             -g {params.input_bins}
#         """

# rule drep_new_ALL_BINS:
#     input: 
#         "data/omics/metagenomes/coassembly/bins/bins_for_drep_ALL_SAMPLES/.bins_linked"
#     output:
#         touch("data/omics/metagenomes/coassembly/bins/.drep_ALL_SAMPLES_done")
#     params:
#         bins_linked = "data/omics/metagenomes/coassembly/bins/bins_for_drep_ALL_SAMPLES/.bins_linked",
#         input_bins = "data/omics/metagenomes/coassembly/bins/bins_for_drep_ALL_SAMPLES/*.fa",
#         genome_info = "data/omics/metagenomes/coassembly/bins/bins_for_drep_ALL_SAMPLES/genome_info.csv", # Will need to be moved to inputs 
#         main_dir = directory("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/"),
#         MAGs= directory("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes")
#     conda: "config/conda_yaml/drep.yaml"
#     resources: cpus=24, mem_mb=150000, time_min=2880
#     shell:
#         """
#         rm -rf {params.main_dir} # Clear any old drep output
        
#         dRep dereplicate \
#             {params.main_dir} \
#             -p {resources.cpus} \
#             --contamination 50 \
#             --completeness 30 \
#             --length 10000 \
#             --genomeInfo {params.genome_info} \
#             -g {params.input_bins}
#         """


# rule mag_coverage_ALL_SAMPLES:
#     input:
#         reads = expand("data/omics/metagenomes/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz", sample = metaG_samples, dir = ["fwd", "rev"]),
#         drep_done = "data/omics/metagenomes/coassembly/bins/.drep_ALL_SAMPLES_done"
#     params:
#         bins_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes"
#     output: 
#         bin_coverage = "data/omics/metagenomes/coassembly/bins/ALL_BINS_coverage.tsv",
#         tmp_dir = temp(directory("data/omics/metagenomes/coassembly/bins/temp_coverm"))
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=24, mem_mb=150000, time_min=2880
#     shell:
#         """
#         # seems to require a temporary directory with lots of space
#         mkdir -p {output.tmp_dir}
#         export TMPDIR={output.tmp_dir}

#         coverm genome \
#             -t {resources.cpus} \
#             -m relative_abundance mean covered_bases variance length rpkm tpm \
#             --min-covered-fraction 0 \
#             --coupled {input.reads} \
#             --genome-fasta-files {params.bins_dir}/*.fa \
#             -o {output.bin_coverage}
#         """

# rule mag_coverage_ALL_SAMPLES_competitive:
#     input:
#         reads = expand("data/omics/metagenomes/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz", sample = metaG_samples, dir = ["fwd", "rev"]),
#         drep_done = "data/omics/metagenomes/coassembly/bins/.drep_ALL_SAMPLES_done"
#     params:
#         bins_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes"
#     output: 
#         bin_coverage = "data/omics/metagenomes/coassembly/bins/ALL_BINS_coverage_competitive.tsv",
#         tmp_dir = temp(directory("data/omics/metagenomes/coassembly/bins/temp_coverm_competitive"))
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=24, mem_mb=150000, time_min=2880
#     shell:
#         """
#         # seems to require a temporary directory with lots of space
#         mkdir -p {output.tmp_dir}
#         export TMPDIR={output.tmp_dir}

#         coverm genome \
#             -t {resources.cpus} \
#             -m relative_abundance mean covered_bases variance length rpkm tpm \
#             --min-covered-fraction 0 \
#             --coupled {input.reads} \
#             --exclude-supplementary \
#             --genome-fasta-files {params.bins_dir}/*.fa \
#             -o {output.bin_coverage}
#         """


# rule GTDB:
#     input:
#         "data/omics/metagenomes/{sample}/bins/bins_for_drep/.bins_linked",
#         "data/reference/GTDBtk/.done_gtdb_refs_downloaded"
#     params:
#         input_bin_dir = "data/omics/metagenomes/{sample}/bins/bins_for_drep",
#         refs = "data/reference/GTDBtk/release207_v2",
#         out_dir = "data/omics/metagenomes/{sample}/bins/GTDB",
#         pplacer_cpus = 1
#     output:
#         done = touch("data/omics/metagenomes/{sample}/bins/.done_GTDB")
#     conda: "config/conda_yaml/gtdbtk.yaml"
#     benchmark: "benchmarks/GTDB/{sample}.txt"
#     log: "logs/GTDB/{sample}.log"
#     resources: cpus=16, mem_mb=100000, time_min=2880
#     shell:
#         """
#         export GTDBTK_DATA_PATH={params.refs}

#         gtdbtk classify_wf \
#             --extension fa \
#             --genome_dir {params.input_bin_dir} \
#             --out_dir {params.out_dir} \
#             --cpus {resources.cpus} \
#             --pplacer_cpus {params.pplacer_cpus}
#         """

# rule run_GTDB:
#     input:
#         expand("data/omics/metagenomes/{sample}/bins/.done_GTDB", sample =metaG_samples)


# rule gunc_GTDB_db_download:
#     output: directory("data/reference/gunc_gtdb")
#     resources: cpus=1, time_min=2880
#     conda: "config/conda_yaml/gunc.yaml"
#     shell:
#         """
#         mkdir -p {output}
#         gunc download_db -db gtdb {output}
#         """

# rule gunc:
#     input:
#         "data/omics/metagenomes/{sample}/bins/bins_for_drep/.bins_linked",
#         ref_file = "data/reference/gunc_gtdb/gunc_db_gtdb95.dmnd",
#         ref_dir = "data/reference/gunc_gtdb"
#     output: 
#         done = touch("data/omics/metagenomes/{sample}/bins/.done_gunc")
#     params: 
#         bin_dir = "data/omics/metagenomes/{sample}/bins/bins_for_drep",
#         out_dir = "data/omics/metagenomes/{sample}/bins/gunc"
#     resources: cpus = 24, mem_mb = 120000, time_min = 2880
#     conda: "config/conda_yaml/gunc.yaml"
#     benchmark: "benchmarks/gunc/{sample}.txt"
#     log: "logs/gunc/{sample}.log"
#     shell:
#         """
#         mkdir -p {params.out_dir}
#         gunc run \
#             --input_dir {params.bin_dir} \
#             -r {input.ref_file} \
#             --threads {resources.cpus} \
#             --temp_dir /tmp \
#             --out_dir {params.out_dir}
#         """
# ####
# # Older binning workflow below
# ####


# # map_reads_to_contigs superseeded by rule map_to_contigs 

# # rule map_reads_to_contigs:
# #     input:
# #         fwd_reads = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
# #         rev_reads = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
# #         contigs = rules.MEC.output.corrected_assembly
# #     output:
# #         unfiltered_bam_dir = temp(directory("data/omics/metagenomes/{sample}/bins/mapped_reads_unfiltered")),
# #         bam = temp("data/omics/metagenomes/{sample}/bins/reads_mapped_to_contigs.bam")
# #     conda: "config/conda_yaml/coverm.yaml"
# #     resources: cpus=24, mem_mb=100000, time_min=2880
# #     shell:
# #         """
# #         coverm make \
# #             -t {resources.cpus} \
# #             -1 {input.fwd_reads} \
# #             -2 {input.rev_reads} \
# #             --reference {input.contigs} \
# #             -o {output.unfiltered_bam_dir}

# #         mv {output.unfiltered_bam_dir}/*.bam {output.bam}
# #         """



# rule dastool:
#     input:
#         contigs = "data/omics/metagenomes/{sample}/bins/work_files/assembly.fa",
#         #bin_folder = rules.make_bins.params.out_dir
#         bins = "data/omics/metagenomes/{sample}/bins/metabat2_bins"
#     params:
#         bin_folder = rules.make_bins.params.out_dir
#     output: 
#         #summary = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_summary.txt",
#         das_folder = directory("data/omics/metagenomes/{sample}/bins/DASTool"),
#         das_bins_folder = directory("data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_bins")
#     conda: "config/conda_yaml/das_tool.yaml"
#     resources: cpus=8, mem_mb=50000, time_min=2880, mem_gb = 50
#     shell:
#         """
#         touch {params.bin_folder}/ran_dastool.touch # for tracking that DAStool ran, even if unsuccessfully

#         #Maxbin2
#         Fasta_to_Contig2Bin.sh \
#         -i {params.bin_folder}/maxbin2_bins \
#         -e fa \
#         > {params.bin_folder}/maxbin.scaffolds2bin.tsv

#         #CONCOT
#         Fasta_to_Contig2Bin.sh \
#         -i {params.bin_folder}/concoct_bins \
#         -e fa \
#         > {params.bin_folder}/concoct.scaffolds2bin.tsv

#         #Metabat2
#         Fasta_to_Contig2Bin.sh \
#         -i {params.bin_folder}/metabat2_bins \
#         -e fa \
#         > {params.bin_folder}/metabat2.scaffolds2bin.tsv

#         DAS_Tool \
#             -i {params.bin_folder}/maxbin.scaffolds2bin.tsv,{params.bin_folder}/concoct.scaffolds2bin.tsv,{params.bin_folder}/metabat2.scaffolds2bin.tsv \
#             -l maxbin,concoct,metabat2 \
#             -c {input.contigs} \
#             -t {resources.cpus} \
#             --write_bins \
#             -o {output.das_folder}/

#         # Rename bins to include sample name, make several downstream analyses easier
#         cd {output.das_bins_folder}
#         for f in *.fa; do mv -v -- "$f" "{wildcards.sample}_$f"; done
#         """


# rule checkm:
#     input: rules.dastool.output.das_bins_folder
#     output:
#         dir = temp(directory("data/omics/metagenomes/{sample}/bins/checkm")),
#         results = "data/omics/metagenomes/{sample}/bins/checkm.txt"
#     conda: "config/conda_yaml/checkm.yaml"
#     resources: cpus=8, mem_mb=80000, time_min=2880, mem_gb = 80
#     shell:
#         """
#         checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {input} {output.dir}
#         """

# rule download_gtdbtk_refs:
#     output:
#         #dir = directory("/home/kiledal/geomicro_home/references/gtdbtk"),
#         #tar = "/home/kiledal/geomicro_home/references/gtdbtk/gtdbtk_data.tar.gz",
#         touch("data/reference/GTDBtk/.done_gtdb_refs_downloaded")
#     params:
#         dir = directory("data/reference/GTDBtk"),
#         tar = "data/reference/GTDBtk/gtdbtk_data.tar.gz"
#     resources: cpus=1, mem_mb=8000, time_min=2880, mem_gb = 8
#     shell:
#         """
#         rm -rf {params.dir}
#         mkdir -p {params.dir}
#         cd {params.dir}
#         wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz
#         #cp ../backup_gtdbtk_data.tar.gz $PWD/gtdbtk_data.tar.gz
#         #tar xzf gtdbtk_data.tar.gz
#         tar xzf gtdbtk_v2_data.tar.gz
#         """


# rule gtdbtk:
#     input:
#         bins = rules.dastool.output.das_bins_folder,
#         refs = "/home/kiledal/geomicro_home/references/gtdbtk/release207"
#     output: directory("data/omics/metagenomes/{sample}/bins/gtdbtk")
#     conda: "config/conda_yaml/gtdbtk.yaml"
#     resources: cpus=1, mem_mb=500000, time_min=2880, mem_gb = 500, partition = "largemem"
#     shell:
#         """
#         GTDBTK_DATA_PATH={input.refs}

#         gtdbtk classify_wf --extension fa --genome_dir {input.bins} --out_dir {output} --cpus {resources.cpus}
#         """

# rule mag_coverage:
#     input:
#         fwd_reads = rules.remove_poly_pairs.output.cleaned_fwd,
#         rev_reads = rules.remove_poly_pairs.output.cleaned_rev,
#         bins = rules.dastool.output.das_bins_folder
#     output: "data/omics/metagenomes/{sample}/bins/coverage.tsv"
#     conda: "config/conda_yaml/coverm_env.yaml"
#     resources: cpus=16, mem_mb=250000, time_min=2880
#     shell:
#         """
#         coverm genome \
#             -t {resources.cpus} \
#             -m relative_abundance mean covered_bases variance length \
#             --min-covered-fraction 0 \
#             -1 {input.fwd_reads} \
#             -2 {input.rev_reads} \
#             --genome-fasta-files {input.bins}/*.fa \
#             -o {output}
#         """

# rule map_to_bins_for_reassembly:
#     input:
#         fwd_reads = rules.remove_poly_pairs.output.cleaned_fwd,
#         rev_reads = rules.remove_poly_pairs.output.cleaned_rev,
#         bin = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_bins/{bin}.fa"
#     output:
#         #unfiltered_bam = "data/omics/metagenomes/{sample}/bins/reassembly/reads/{bin}_prefilt.bam",
#         unfiltered_bam_dir = temp(directory("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}__bam_unfiltered")),
#         filtered_bam = temp("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}.bam"),
#         #filtered_bam_dir = directory("data/omics/metagenomes/{sample}/bins/reassembly/reads/{bin}__bam_filtered"),
#         fwd_reads = temp("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}_R1.fastq.gz"),
#         rev_reads = temp("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}_R2.fastq.gz")
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=24, mem_mb=100000, time_min=2880
#     shell:
#         """
#         coverm make \
#             -t {resources.cpus} \
#             -1 {input.fwd_reads} \
#             -2 {input.rev_reads} \
#             --reference {input.bin} \
#             -o {output.unfiltered_bam_dir}

#         coverm filter \
#             -b {output.unfiltered_bam_dir}/*.bam \
#             -o {output.filtered_bam} \
#             --min-read-percent-identity 0.98 \
#             --threads {resources.cpus}

#         samtools fastq -@ {resources.cpus} \
#             {output.filtered_bam} \
#             -1 {output.fwd_reads} \
#             -2 {output.rev_reads} \
#             -0 /dev/null -s /dev/null -n
#         """

# rule reassemble_bin:
#     input:
#         fwd_reads = rules.map_to_bins_for_reassembly.output.fwd_reads,
#         rev_reads = rules.map_to_bins_for_reassembly.output.rev_reads,
#         bin = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_bins/{bin}.fa"
#     output:
#         #assembly_dir = directory("data/omics/metagenomes/{sample}/bins/reassembly/{bin}"),
#         contigs = "data/omics/metagenomes/{sample}/bins/reassembly/{bin}_reassembled_contigs.fasta"
#     params: 
#         assembly_dir = directory("data/omics/metagenomes/{sample}/bins/reassembly/{bin}")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/bin_reassembly/{sample}/{bin}.log"
#     benchmark: "logs/bin_reassembly/{sample}/{bin}.log"
#     resources: cpus = 16, time_min=20160, mem_mb = 16000
#         #mem_mb = lambda wildcards, attempt: attempt * 100000
#     shell:
#         """
#         spades.py \
#             -t {resources.cpus} \
#             -m $(echo "scale=-1; ({resources.mem_mb}/1000)/1" | bc) \
#             --careful \
# 		    --untrusted-contigs {input.bin} \
# 		    -1 {input.fwd_reads} \
# 		    -2 {input.rev_reads} \
# 		    -o {params.assembly_dir} > {log}

#         mv {params.assembly_dir}/contigs.fasta {output.contigs}
#         rm -r {params.assembly_dir}
#         """


# rule gather_bin:
#     input: expand("data/omics/metagenomes/{sample}/bins/ran_dastool.touch",sample=metaG_samples)
#     output: directory("data/omics/metagenome_bins")
#     shell:
#         """
#         mkdir -p {output}
        
#         ln data/omics/metagenomes/*/bins/DASTool/_DASTool_bins/*.fa {output}
#         """

# # rule gtdbtk_all:
# #     input:
# #         bins = rules.gather_bin.output,
# #         refs = "/geomicro/data2/kiledal/references/gtdbtk/release202"
# #     output: directory("data/gtdbtk")
# #     conda: "code/gtdbtk.yaml"
# #     resources: cpus=32, mem_mb=250000, time_min=2880, mem_gb = 250
# #     shell:
# #         """
# #         GTDBTK_DATA_PATH={input.refs}

# #         gtdbtk classify_wf --extension fa --genome_dir {input.bins} --out_dir {output} --cpus {resources.cpus} --pplacer_cpus 1
# #         """

# checkpoint drep:
#     #input: 
#     #    bins = rules.gather_bin.output
#     output:
#         main_dir = directory("data/omics/metagenomes/metagenome_bins/derep"),
#         MAGs= directory("data/omics/metagenomes/metagenome_bins/derep/dereplicated_genomes")
#     conda: "config/conda_yaml/drep.yaml"
#     resources: cpus=8, mem_mb=250000, time_min=2880, partition = "largemem"
#     shell:
#         """
#         dRep dereplicate \
#             {output.main_dir} \
#             -p {resources.cpus} \
#             -g data/omics/metagenomes/metagenome_bins/*.fa
#         """

# rule drep_mag_coverage:
#     input:
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
#         bins = "data/omics/metagenomes/metagenome_bins/derep/dereplicated_genomes"
#     output: "data/omics/metagenomes/{sample}/bins/drep_bins_coverage.tsv"
#     conda: "config/conda_yaml/coverm.yaml"
#     resources: cpus=36, mem_mb=170000, time_min=2880
#     shell:
#         """
#         coverm genome \
#             -t {resources.cpus} \
#             -m relative_abundance mean covered_bases variance length \
#             --min-covered-fraction 0 \
#             -1 {input.fwd_reads} \
#             -2 {input.rev_reads} \
#             --genome-fasta-files {input.bins}/*.fa \
#             -o {output}
#         """

# rule run_drep_bin_coverage:
#     input:
#         expand("data/omics/metagenomes/{sample}/bins/drep_bins_coverage.tsv", sample = metaG_samples)


# rule get_traitar_pfamDB:
#     output: directory("data/reference/traitar_pfamDB")
#     container: "library://a_gihawi/traitar3/traitar3"
#     resources: cpus = 1, mem_mb = 20000
#     shell:
#         """
#         traitar pfam {output}
#         """

# rule traitar:
#     input:
#         #bin_genes = expand("data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/metagenome_bins/{BIN,[^/]+}.fa").BIN),
#         pfam_db = "data/reference/traitar_pfamDB",
#         #bin_dir = "data/omics/metagenomes/metagenome_bins",
#         sample_file = "data/omics/metagenomes/metagenome_bins/traitar_sample_list.tsv"
#     params:
#         bin_dir = "data/omics/metagenomes/metagenome_bins/derep/dereplicated_genomes",
#         gene_dir = "data/omics/metagenomes/metagenome_bins/prodigal",
#         out_dir = "data/omics/metagenomes/metagenome_bins/traitar"
#     output: 
#         #directory("data/omics/metagenomes/metagenome_bins/traitar")
#         touch("data/omics/metagenomes/metagenome_bins/traitar/.done")
#     container: "library://a_gihawi/traitar3/traitar3"
#     resources: cpus = 1, mem_mb = 170000, time_min=2880 #, partition = "largemem"
#     shell:
#         """
#         #mkdir -p {params.out_dir}

#         traitar phenotype --overwrite -c {resources.cpus} /db {params.gene_dir} {input.sample_file} from_genes {params.out_dir}
#         """

# rule fegenie:
#     input:
#         #bin_genes = expand("data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/metagenome_bins/{BIN,[^/]+}.fa").BIN),
#         bin_dir = "data/omics/metagenomes/metagenome_bins",
#     output: directory("data/omics/metagenomes/metagenome_bins/fegenie")
#     conda: "config/conda_yaml/fegenie.yaml"
#     resources: cpus = 16, mem_mb = 50000, time_min=2880
#     shell:
#         """
#         FeGenie.py -bin_dir data/omics/metagenomes/metagenome_bins/prodigal/ -bin_ext faa -t {resources.cpus} -out {output} --orfs
#         """


# rule prodigal_mags:
#     input:
#         bin = "data/omics/metagenomes/metagenome_bins/{bin}.fa"
#     output: 
#         proteins = "data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa",
#         genes = "data/omics/metagenomes/metagenome_bins/prodigal/{bin}.gff"
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus = 1, mem_mb = 10000
#     shell:
#         """
#         prodigal -p meta -i {input.bin} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
#         """

# rule run_prodigal_mags:
#     input: expand("data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/metagenome_bins/{BIN,[^/]+}.fa").BIN)

# rule prodigal_mags_DREP:
#     input:
#         bin = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{bin}.fa"
#     output: 
#         proteins = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa",
#         genes = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.gff"
#     conda: "config/conda_yaml/main.yaml"
#     resources: cpus = 1, mem_mb = 10000
#     shell:
#         """
#         prodigal -p meta -i {input.bin} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
#         """

# rule run_prodigal_mags_DREP:
#     input: expand("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{BIN,[^/]+}.fa").BIN)


# rule traitar_coassembly:
#     input:
#         #bin_genes = expand("data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/metagenome_bins/{BIN,[^/]+}.fa").BIN),
#         pfam_db = "data/reference/traitar_pfamDB",
#         #bin_dir = "data/omics/metagenomes/metagenome_bins",
#         sample_file = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar_sample_list.tsv"
#     params:
#         bin_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes",
#         gene_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal",
#         out_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar"
#     output: 
#         #directory("data/omics/metagenomes/metagenome_bins/traitar")
#         touch("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar/.done")
#     benchmark: "logs/traitar_coassembly/benchmark.txt"
#     container: "library://a_gihawi/traitar3/traitar3"
#     resources: cpus = 16, mem_mb = 170000, time_min=20000 #, partition = "largemem"
#     shell:
#         """
#         #mkdir -p {params.out_dir}

#         pwd
#         cd ~/scratch_gdick1/GVHD/
#         pwd

#         traitar phenotype --overwrite -c {resources.cpus} /db {params.gene_dir} {input.sample_file} from_genes {params.out_dir}
#         """

# rule fegenie_coassembly:
#     input:
#         #bin_genes = expand("data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/metagenome_bins/{BIN,[^/]+}.fa").BIN),
#         bin_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes"
#     params:
#         prodigal_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal"
#     output: directory("data/omics/metagenomes/coassembly/bins/fegenie")
#     conda: "config/conda_yaml/fegenie.yaml"
#     resources: cpus = 16, mem_mb = 50000, time_min=2880
#     shell:
#         """
#         FeGenie.py -bin_dir {params.prodigal_dir}/ -bin_ext faa -t {resources.cpus} -out {output} --orfs
#         """


# rule kofam_scan_bins:
#     input:
#         genes = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa",
#         profile = "data/reference/kegg/kofamscan/profiles",
#         ko_list = "data/reference/kegg/kofamscan/ko_list"
#     output:
#         ko_annot = "data/omics/metagenomes/coassembly/bins/kofamscan/bins/{bin}_kofam_results.txt"
#     conda: "config/conda_yaml/kofamscan.yaml"
#     #shadow: "shallow"
#     benchmark: "benchmarks/kofamscan/{bin}.txt"
#     log: "logs/kofamscan/{bin}.log"
#     resources: cpus=12, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 50000
#     shell:
#         """
#         exec_annotation \
#             -o {output.ko_annot} \
#             --format=detail-tsv \
#             --cpu={resources.cpus}  \
#             --profile {input.profile} \
#             --tmp-dir=/tmp/{wildcards.bin}_kofamscan \
#             --ko-list {input.ko_list} {input.genes}
#         """

# rule minpath_bins:
#     input: 
#         minpath = "code/MinPath/MinPath.py",
#         KOs = "data/omics/metagenomes/coassembly/bins/minpath/input/{bin}.txt"
#     output: 
#         report = "data/omics/metagenomes/coassembly/bins/minpath/output/{bin}.txt",
#         detailed = "data/omics/metagenomes/coassembly/bins/minpath/output/{bin}__detailed.txt"
#     params:
#         mps = "data/omics/metagenomes/coassembly/bins/minpath/output/{bin}.mps"
#     resources: cpus = 1, mem_mb = 16000, time = 2000
#     shell:
#         """
#         python {input.minpath} -ko {input.KOs} -report {output.report} -details {output.detailed} -mps {params.mps}
#         """

# rule minpath_bins_ec:
#     input: 
#         minpath = "code/MinPath/MinPath.py",
#         ec = "data/omics/metagenomes/coassembly/bins/minpath/input_ec/{bin}.txt"
#     output: 
#         report = "data/omics/metagenomes/coassembly/bins/minpath/output_ec/{bin}.txt",
#         detailed = "data/omics/metagenomes/coassembly/bins/minpath/output_ec/{bin}__detailed.txt"
#     params:
#         mps = "data/omics/metagenomes/coassembly/bins/minpath/output_ec/{bin}.mps"
#     resources: cpus = 1, mem_mb = 16000, time = 2000
#     shell:
#         """
#         python {input.minpath} -ec {input.ec} -report {output.report} -details {output.detailed} -mps {params.mps}
#         """

# rule run_minpath_bins_ec:
#     input:
#         expand("data/omics/metagenomes/coassembly/bins/minpath/output_ec/{bin}.txt", bin = glob_wildcards("data/omics/metagenomes/coassembly/bins/minpath/input_ec/{BIN,[^/]+}.txt").BIN)


# rule kofam_scan_assemblies:
#     input:
#         genes = "data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa",
#         profile = "data/reference/kegg/kofamscan/profiles",
#         ko_list = "data/reference/kegg/kofamscan/ko_list"
#     output:
#         ko_annot = "data/omics/metagenomes/{sample}/proteins/{sample}_kofam_results.txt"
#     conda: "config/conda_yaml/kofamscan.yaml"
#     #shadow: "shallow"
#     benchmark: "benchmarks/kofamscan_assembly/{sample}.txt"
#     log: "logs/kofamscan/{sample}.log"
#     resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 120000
#     shell:
#         """
#         exec_annotation \
#             -o {output.ko_annot} \
#             --format=detail-tsv \
#             --cpu={resources.cpus}  \
#             --profile {input.profile} \
#             --tmp-dir=/tmp/{wildcards.sample}_kofamscan \
#             --ko-list {input.ko_list} {input.genes}
#         """

# rule minpath_assembly:
#     input: 
#         minpath = "code/MinPath/MinPath.py",
#         KOs = "data/omics/metagenomes/{sample}/minpath/KOs.txt"
#     output: 
#         report = "data/omics/metagenomes/{sample}/minpath/{sample}__report.txt",
#         detailed = "data/omics/metagenomes/{sample}/minpath/{sample}__detailed.txt"
#     params:
#         mps = "data/omics/metagenomes/{sample}/minpath/{sample}.mps"
#     resources: cpus = 1, mem_mb = 16000, time = 20000
#     shell:
#         """
#         python {input.minpath} -ko {input.KOs} -report {output.report} -details {output.detailed} -mps {params.mps}
#         """


# rule run_kofam_scan:
#     input:
#         expand("data/omics/metagenomes/coassembly/bins/kofamscan/bins/{bin}_kofam_results.txt", bin = glob_wildcards("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{BIN,[^/]+}.fa").BIN),
#         expand("data/omics/metagenomes/coassembly/bins/minpath/output/{bin}.txt", bin = glob_wildcards("data/omics/metagenomes/coassembly/bins/minpath/input/{BIN,[^/]+}.txt").BIN)

# rule run_kofam_scan_assemblies:
#     input:
#         expand("data/omics/metagenomes/{sample}/proteins/{sample}_kofam_results.txt", sample = metaG_samples),
#         expand("data/omics/metagenomes/{sample}/minpath/{sample}__report.txt", sample = glob_wildcards("data/omics/metagenomes/{samp}/minpath/KOs.txt").samp)

# # rule humann_db_download:
# #     output:
# #         "data/reference/humann/all_genes_annot.1.bt2l",
# #         "data/reference/humann/all_genes_annot.2.bt2l",
# #         "data/reference/humann/all_genes_annot.3.bt2l",
# #         "data/reference/humann/all_genes_annot.4.bt2l",
# #         "data/reference/humann/all_genes_annot.rev.1.bt2l",
# #         "data/reference/humann/all_genes_annot.rev.2.bt2l",
# #         "data/reference/humann/genome_reps_filt_annot.faa.gz",
# #         "data/reference/humann/genome_reps_filt_annot.fna.gz",
# #         "data/reference/humann/genome_reps_filt_annot.tsv.gz",
# #         "data/reference/humann/protein_database/uniref90_201901b.dmnd"
# #     resources: time_min=1440,
# #     shell:
# #         """
# #         cd data/reference/humann/
# #         wget -nv -r -nH --cut-dirs=6 -nc ftp://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release95/humann3/uniref90/

# #         #Some naming issues with the diamond database, these may get cleared up over time. But for now...
# #         mv -u data/reference/humann/protein_database/uniref90_201901.dmnd data/reference/humann/protein_database/uniref90_201901b.dmnd
# #         mv -u data/reference/humann/protein_database/uniref90_201901.md5 data/reference/humann/protein_database/uniref90_201901b.md5
# #         """


# rule install_metaphlan_db:
#     output: directory("data/reference/metaphlan")
#     conda:
#         "config/conda_yaml/metaphlan.yaml"
#     resources:
#         cpus = 24, mem_mb = 150000, time_min=10000
#     shell:
#         """
#         metaphlan --install --bowtie2db {output}
#         """

# rule metaphlan:
#     input:
#         fwd_reads = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         rev_reads = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
#         refs = "data/reference/metaphlan"
#     output:
#         abund = "data/omics/metagenomes/{sample}/metaphlan/{sample}_abund.txt",
#         sam = "data/omics/metagenomes/{sample}/metaphlan/{sample}.sam.bz2",
#         bt2 = "data/omics/metagenomes/{sample}/metaphlan/{sample}.bowtie2.bz2"
#     params:
#     conda:
#         "config/conda_yaml/metaphlan.yaml"
#     log:
#         "logs/metaphlan/{sample}.log"
#     benchmark:
#         "benchmarks/metaphlan/{sample}.txt"
#     resources:
#         cpus = 24, mem_mb = 150000, time_min=10000
#     shell:
#         """
#         metaphlan \
#             {input.fwd_reads},{input.rev_reads} \
#             -o {output.abund} \
#             -s {output.sam} \
#             --bowtie2out {output.bt2} \
#             --input_type fastq \
#             --unclassified_estimation \
#             --add_viruses \
#             --bowtie2db {input.refs} \
#             --nproc {resources.cpus} | tee {log}
#         """

# rule run_metaphlan:
#     input: expand("data/omics/metagenomes/{sample}/metaphlan/{sample}_abund.txt", sample = metaG_samples)

# rule humann:
#     input:
#         # Sample data #
#         f_seq = rules.remove_poly_pairs.output.cleaned_fwd,
#         r_seq = rules.remove_poly_pairs.output.cleaned_rev,
#         bracken_mpa = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_brackenMpa.txt",
#         # Reference database #
#         NUC_DB = "data/reference/humann/genome_reps_filt_annot.fna.gz",
#         PROT_DB = "data/reference/humann/protein_database/uniref90_201901b.dmnd",
#         NUC_fol = "data/reference/humann/",
#         PROT_fol = "data/reference/humann/protein_database/"
#     output:
#         humann_output = directory("data/omics/metagenomes/{sample}/humann"),
#         concat_unzipped_reads = temp("data/omics/metagenomes/{sample}/reads/for_humann.fastq")
#     params:
#         mem_use = "maximum"
#     conda: "config/conda_yaml/humann.yaml"
#     log: "logs/humann/{sample}.log"
#     benchmark: "benchmarks/humann/{sample}.txt"
#     #resources: cpus=36, mem_mb=175000, time_min=20160
#     resources: cpus=36, time_min=20160, partition = "largemem", mem_mb = lambda wildcards, attempt: attempt * 250000
#     shell:
#         """
#         #Humann needs non-compressed fastqs, and forward and reverse files should be concatenated
#         gunzip -c {input.f_seq} > {output.concat_unzipped_reads}
#         gunzip -c {input.r_seq} >> {output.concat_unzipped_reads}

#         echo "Combined and uncompressed fastq files.\n"

#         #proj_dir=$PWD

#         humann3 --bypass-nucleotide-index \
#             --threads {resources.cpus} \
#             --memory-use {params.mem_use} \
#             --nucleotide-database {input.NUC_fol} \
#             --protein-database {input.PROT_fol} \
#             --taxonomic-profile {input.bracken_mpa} \
#             --input {output.concat_unzipped_reads} \
#             --output-basename {wildcards.sample}_humann \
#             --output {output.humann_output}/ > {log}
#         """

# rule run_humann: 
#     input: expand("data/omics/metagenomes/{sample}/humann", sample=metaG_samples)

# rule humann_fastp:
#     input:
#         # Sample data #
#         f_seq = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#         r_seq = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
#         bracken_mpa = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_brackenMpa.txt",
#         # Reference database #
#         NUC_DB = "data/reference/humann/genome_reps_filt_annot.fna.gz",
#         PROT_DB = "data/reference/humann/protein_database/uniref90_201901b.dmnd",
#         NUC_fol = "data/reference/humann/",
#         PROT_fol = "data/reference/humann/protein_database/"
#     output:
#         humann_output = directory("data/omics/metagenomes/{sample}/humann_fastp"),
#         concat_unzipped_reads = temp("data/omics/metagenomes/{sample}/reads/for_humann_fastp.fastq")
#     params:
#         mem_use = "maximum"
#     conda: "config/conda_yaml/humann.yaml"
#     log: "logs/humann/{sample}_fastp.log"
#     benchmark: "benchmarks/humann/{sample}_fastp.txt"
#     #resources: cpus=36, mem_mb=175000, time_min=20160
#     #resources: cpus=36, time_min=20160, partition = "largemem", mem_mb = lambda wildcards, attempt: attempt * 250000
#     resources: cpus=36, time_min=20160, mem_mb = 175000
#     shell:
#         """
#         #Humann needs non-compressed fastqs, and forward and reverse files should be concatenated
#         gunzip -c {input.f_seq} > {output.concat_unzipped_reads}
#         gunzip -c {input.r_seq} >> {output.concat_unzipped_reads}

#         echo "Combined and uncompressed fastq files.\n"

#         #proj_dir=$PWD

#         humann3 --bypass-nucleotide-index \
#             --threads {resources.cpus} \
#             --memory-use {params.mem_use} \
#             --nucleotide-database {input.NUC_fol} \
#             --protein-database {input.PROT_fol} \
#             --taxonomic-profile {input.bracken_mpa} \
#             --input {output.concat_unzipped_reads} \
#             --output-basename {wildcards.sample}_humann \
#             --output {output.humann_output}/ > {log}

#         rm -r {output.humann_output}/{wildcards.sample}_humann_humann_temp
#         """

# rule run_humann_fastp: 
#     input: expand("data/omics/metagenomes/{sample}/humann_fastp", sample=metaG_samples)



# rule sourmash_gather:
#     input:
#        f_seq = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
#        r_seq = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
#        refDB = "data/reference/sourmash/gtdb-rs207.dna.k31.zip",
#        taxDB = "data/reference/sourmash/gtdb-rs207.taxonomy.sqldb"
#     output:
#        sig = "data/omics/metagenomes/{sample}/sourmash/{sample}.sig",
#        reps = "data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.csv",
#        tax = "data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.with-lineages.csv"
#     conda: "config/conda_yaml/sourmash.yaml"
#     log: "logs/sourmash/{sample}.log"
#     benchmark: "benchmarks/sourmash/{sample}.txt"
#     resources: cpus=2, time_min=4320, mem_mb = 20000
#     shell:
#         """
#         sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge {wildcards.sample} -o {output.sig} {input.f_seq} {input.r_seq} 2>&1 | tee {log}

#         sourmash gather {output.sig} {input.refDB} -o {output.reps} 2>&1 | tee -a {log}

#         sourmash tax annotate -g {output.reps} -t {input.taxDB} 2>&1 | tee -a {log}

#         mv $(basename {output.tax}) data/omics/metagenomes/{wildcards.sample}/sourmash/
#         """


# rule run_sourmash: 
#     input: expand("data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.csv", sample=metaG_samples)


# rule setup_DRAM:
#     output:
#         db = directory("data/reference/DRAM")
#     conda: "config/conda_yaml/DRAM.yaml"
#     resources: 
#         cpus = 24, partition = "largemem", mem_mb=600000, time_min=2000
#     shell:
#         """
#         mkdir -p {output.db}
#         DRAM-setup.py prepare_databases --output_dir {output.db}
#         """

# rule DRAM_annotate:
#     input: "data/omics/metagenomes/coassembly/bins/.drep_ALL_SAMPLES_done"
#     output:
#         dram_annotation = directory("data/omics/metagenomes/coassembly/bins/DRAM")
#     params:
#         bin_dir= "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes",
#         bin_extension = ".fa"
#     conda: "config/conda_yaml/DRAM.yaml"
#     log: "logs/DRAM_annotate.log"
#     benchmark: "benchmarks/DRAM_annotate.txt"
#     resources: cpus = 24, mem_mb = 500000, time_min=20000, partition = "largemem"
#     shell:
#         """
#         DRAM.py annotate \
#             -i '{params.bin_dir}/*{params.bin_extension}' \
#             -o {output.dram_annotation} \
#             --min_contig_size 1000 \
#             --threads {resources.cpus}
#         """

# rule DRAM_annotate_indiv:
#     input: 
#         genome = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{genome}.fa"
#     output:
#         dram_annotation = directory("data/omics/metagenomes/coassembly/bins/DRAM/indiv_bin_annot/{genome}")
#     conda: "config/conda_yaml/DRAM.yaml"
#     log: "logs/dram_annotate_indiv/{genome}.log"
#     benchmark: "benchmarks/dram_annotate_indiv/{genome}.txt"
#     resources: cpus = 16, mem_mb = 100000, time_min=20000
#     shell:
#         """
#         DRAM.py annotate \
#             -i {input.genome} \
#             -o {output.dram_annotation} \
#             --min_contig_size 1000 \
#             --threads {resources.cpus}
#         """

# rule DRAM_distill:
#     input: 
#         concat_dram_annotation = "data/omics/metagenomes/coassembly/bins/DRAM/indiv_bin_annot_concat/annotations.tsv",
#         concat_dram_trna = "data/omics/metagenomes/coassembly/bins/DRAM/indiv_bin_annot_concat/tRNAs.tsv",
#         concat_dram_rrna = "data/omics/metagenomes/coassembly/bins/DRAM/indiv_bin_annot_concat/rRNAs.tsv",
#     output:
#         distillate = directory("data/omics/metagenomes/coassembly/bins/DRAM/indiv_bin_annot_concat/genome_summaries")
#     conda: "config/conda_yaml/DRAM.yaml"
#     log: "logs/DRAM_distill.log"
#     benchmark: "benchmarks/dram_distill.txt"
#     resources: cpus = 16, mem_mb = 100000, time_min=20000
#     shell:
#         """
#         DRAM.py distill \
#             -i {input.concat_dram_annotation} \
#             -o {output.distillate} \
#             --trna_path {input.concat_dram_trna} \
#             --rrna_path {input.concat_dram_rrna}
#         """


# rule run_DRAM_indiv:
#     input:
#         expand("data/omics/metagenomes/coassembly/bins/DRAM/indiv_bin_annot/{genome}", genome = glob_wildcards("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{genome_name}.fa").genome_name)


# rule bakta:
#     input:
#         genome = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{genome}.fa"
#     output:
#         dir = directory("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/bakta/{genome}")
#     params:
#         db = "data/reference/bakta/db"
#     conda: "config/conda_yaml/bakta.yaml"
#     log: "logs/bakta/{genome}.tsv"
#     benchmark: "benchmarks/bakta/{genome}.tsv"
#     resources: cpus=24, mem_mb=48000, time_min=10000, 
#     shell:
#         """
#         bakta --db {params.db} \
#             --output {output.dir} \
#             --threads {resources.cpus} \
#             {input.genome} | tee {log}
#         """

# import pandas as pd
# bakta_df = pd.read_table('data/omics/metagenomes/coassembly/bins/drep_bakta_dirs.tsv').set_index("bakta_dir", drop=False)
# bakta_genome_names = list(bakta_df['bakta_dir'])

# rule run_bakta:
#     input:
#         bakta_genome_names



# rule pathway_tools:
#     output: touch("data/.pathway_tools_run")
#     conda: "pathwaytools"
#     params:
#         annot_dir="data/omics/metagenomes/coassembly/bins/pathway_tools_inputs",
#         output_dir = "data/omics/metagenomes/coassembly/bins/pathway_tools",
#         tax_file = "data/omics/metagenomes/coassembly/bins/pathway_tools_tax.tsv"
#     resources: cpus=32, mem_mb=90000, time_min=20000
#     shell:
#         """
#         export PATH="$HOME/GVHD/code/pathwaytools/pathway-tools:$PATH"

#         mpwt \
#             -f={params.annot_dir} \
#             -o={params.output_dir} \
#             --patho \
#             --hf \
#             --op \
#             --tp \
#             --nc \
#             --flat \
#             --mc \
#             --cpu={resources.cpus} \
#             --log=logs/pathwaytools \
#             --taxon-file {params.tax_file}
#         """


# rule g_to_tree:
#     output:
#         outdir = directory("data/omics/metagenomes/coassembly/bins/gtotree_dreped")
#     params:
#         genomes_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/",
#         genome_ext = ".fa",
#         genome_list = "/tmp/genomes_for_gtotree.txt"
#     conda: "config/conda_yaml/g_to_tree.yaml"
#     log: "logs/g_to_tree.log"
#     benchmark: "benchmarks/g_to_tree.txt"
#     resources: cpus = 12, mem_mb = 50000, time_min=20000
#     shell: 
#         """
#         ls {params.genomes_dir}/*{params.genome_ext} > {params.genome_list}
        
#         GToTree \
#         -f {params.genome_list} \
#         -o {output.outdir} \
#         -H Bacteria_and_Archaea.hmm \
#         -j {resources.cpus} \
#         -n 1
#        """

rule reads_unirefLCA_mmseqs:
    input:
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        report = "data/omics/{sample_type}/{sample}/uniref_readmapping/{sample}_report"
    conda:  "config/conda_yaml/mmseqs.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/uniref_readmapping/{sample}",
        # tmp_dir = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/tmp",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        tmp_dir = "/tmp/kiledal/mmseqs2/{sample}",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/reads_unirefLCA_mmseqs/{sample_type}-{sample}.txt"
    log: "logs/reads_unirefLCA_mmseqs/{sample_type}-{sample}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        mem_mb = 160000, cpus=32, time_min=20000
    shell:
        """
        mkdir -p {params.tmp_dir}
        mkdir -p $(dirname {params.out_prefix})

        touch {params.out_prefix}_test.touch

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h {params.tmp_dir} 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        # Previous copied reads to temp space, but no longer think this is nescessary. Should be deleted.
        #cp {input.fwd_reads} {params.tmp_fwd_reads}
        #cp {input.rev_reads} {params.tmp_rev_reads}

        #mmseqs touchdb {params.unirefDB} # loads  the database into memory, only use in large memory systems / nodes
        
        mmseqs \
            easy-taxonomy \
            {input.fwd_reads} {input.rev_reads} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            {params.tmp_dir} \
            --lca-mode 3 \
            --orf-filter 1 \
            --orf-filter-s 3.5 \
            -s 4 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit 100G \
            2>&1 | tee -a {log}

            #{params.tmp_fwd_reads} {params.tmp_rev_reads} \
            # --db-load-mode 2 # this loads the database into memory, was an attempt to speed up on Great Lakes but limited memory & scratch space were limiting factors

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h {params.tmp_dir} 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r {params.tmp_fwd_reads} {params.tmp_rev_reads} {params.tmp_dir} 2>&1 | tee -a {log}
        rm -r {params.tmp_dir} 2>&1 | tee -a {log}
        printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """


rule add_lineage_to_unirefLCAtax:
    input:
        db = ancient("/home/kiledal/geomicro_home/GLAMR/data/reference/ncbi_tax"),
        report = "data/omics/metagenomes/{sample}/uniref_readmapping/{sample}_report"
    output: 
        inspect_w_lineage = "data/omics/metagenomes/{sample}/uniref_readmapping/{sample}_report_w_full_lineage",
        inspect_w_7_lev_lineage = "data/omics/metagenomes/{sample}/uniref_readmapping/{sample}_report_w_standardized_lineage"
    conda: "config/conda_yaml/taxonkit.yaml"
    resources: cpus=1, mem_mb=10000, time_min=5440, mem_gb = 10
    shell:
        """
        taxonkit lineage \
            {input.report} \
            --data-dir {input.db} \
            -i 5 \
            -o {output.inspect_w_lineage} 2>&1 | tee {log}

        taxonkit reformat \
            {output.inspect_w_lineage} \
            --data-dir {input.db} \
            -i 7 \
            -P \
            -o {output.inspect_w_7_lev_lineage} 2>&1 | tee -a {log}
        """


rule run_mmseqsLCA_GL:
    input: 
        expand("data/omics/metagenomes/{sample}/uniref_readmapping/{sample}_report", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/uniref_readmapping/{sample}_report_w_full_lineage", sample = metaG_samples)


rule contig_abund:
    input:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit/final.contigs.fa",
        f_seq = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output: 
        coverage_full = "data/omics/{sample_type}/{sample}/{sample}_contig_abund.tsv"
    params:
        tmpdir = "tmp/coverm_contig_abund/{sample}"
    conda: "config/conda_yaml/coverm.yaml"
    log: "logs/contig_abund/{sample}-{sample_type}.log"
    benchmark: "benchmarks/contig_abund/{sample}-{sample_type}.txt"
    resources: cpus=24, mem_mb=120000, time_min=360 # standard assemblies
    #resources: cpus=24, mem_mb=1000000, time_min=2880, partition = "largemem" # coassembly
    shell:
        """
        mkdir -p {params.tmpdir}
        TMPDIR={params.tmpdir}

        # Link reads w/ naming convention prefered by coverM
        fwd_reads=$(dirname {input.f_seq})/{wildcards.sample}_R1.fq.gz
        rev_reads=$(dirname {input.r_seq})/{wildcards.sample}_R2.fq.gz
        ln {input.f_seq} $fwd_reads
        ln {input.r_seq} $rev_reads

        coverm contig \
            -c $fwd_reads $rev_reads \
            -r {input.contigs} \
            -t {resources.cpus} \
            -m mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
            --output-format sparse \
            --output-file {output.coverage_full} 2>&1 | tee {log}

        rm -r {params.tmpdir} $fwd_reads $rev_reads 
        """


rule contig_unirefLCA_mmseqs:
    input:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit/final.contigs.fa"
    output:
        report = "data/omics/{sample_type}/{sample}/{sample}_contig_report"
    conda:  "config/conda_yaml/mmseqs.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/{sample}_contig",
        # tmp_dir = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/tmp",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        tmp_dir = "/tmp/kiledal/mmseqs2/{sample}_contigs",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/contig_unirefLCA_mmseqs/{sample_type}-{sample}.txt"
    log: "logs/contig_unirefLCA_mmseqs/{sample_type}-{sample}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        mem_mb = 160000, cpus=32, time_min=1440
    shell:
        """
        mkdir -p {params.tmp_dir}

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h {params.tmp_dir} 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #mmseqs touchdb {params.unirefDB} # loads the database into memory, only use in large memory systems / nodes

        mmseqs \
            easy-taxonomy \
            {input.contigs} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            {params.tmp_dir} \
            --lca-mode 3 \
            --orf-filter 1 \
            --orf-filter-s 3.5 \
            -s 4 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit $(echo "scale=0;({resources.mem_mb}*0.8)/1024" | bc)G \
            --db-load-mode 1 \
            2>&1 | tee -a {log}

            # --db-load-mode 2 # this loads the database into memory, was an attempt to speed up on Great Lakes but limited memory & scratch space were limiting factors

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h {params.tmp_dir} 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r {params.tmp_fwd_reads} {params.tmp_rev_reads} {params.tmp_dir} 2>&1 | tee -a {log}
        rm -r {params.tmp_dir} 2>&1 | tee -a {log}
        printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """

rule tax_abund_summary_from_contigs:
    input: 
        mmseqs_report = "data/omics/{sample_type}/{sample}/{sample}_contig_report",
        script = "code/tax_abund_from_contigs.R",
        contig_abund = "data/omics/{sample_type}/{sample}/{sample}_contig_abund.tsv"
        #assembly_done = "data/omics/{sample_type}/{sample}/assembly/megahit/.done"
    output:
        abund_summary = "data/omics/{sample_type}/{sample}/{sample}_lca_abund_summarized.tsv"
    params:
        lca = "data/omics/{sample_type}/{sample}/{sample}_contig_lca.tsv",
        taxonkit_path = "~/miniconda3/envs/taxonkit/bin/taxonkit",
        taxdump = "data/reference/ncbi_tax"
    benchmark: "benchmarks/tax_abund_summary_from_contigs/{sample_type}-{sample}.txt"
    container: "docker://eandersk/r_microbiome"
    resources: cpus = 24, time_min=1000, mem_mb = 100000
    shell:
        """
        pwd
        
        ./{input.script} \
            -l {params.lca} \
            -r {input.contig_abund} \
            -o {output.abund_summary} \
            -c {resources.cpus} \
            -t {params.taxonkit_path} \
            -d {params.taxdump}
        """


rule run_mmseqsLCA_contig_largemem:
    input: 
        expand("data/omics/metagenomes/{sample}/{sample}_contig_report", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/{sample}_contig_abund.tsv", sample = metaG_samples),
        expand("data/omics/{sample_type}/{sample}/{sample}_lca_abund_summarized.tsv", sample = metaG_samples, sample_type = "metagenomes")