cd /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/fasta_inputs/Share_Heidelberg/Oct2024/
for sample in `ls`; do
        mkdir -p /hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/ubam2counts/ubam2counts_HD_Oct2024/$sample
        output_file=/hpc/pmc_kool/fvalzano/wdl_pipeline_v12.1.0/wdl/inputs/ubam2counts/ubam2counts_HD_Oct2024/$sample/rna_fusions_germline_snv_inputs_no_molgenis_fv.json
        # Create json content, this is the standard form, without the sample specific info that you will have to enter manually
        json_content='{
  "rnafusion.ApplyBQSR_RNA.emit_oq_flag": "true",
  "rnafusion.ApplyBQSR_RNA.sequence_group_interval": null,
  "rnafusion.BaseRecalibrator_RNA.sequence_group_interval": null,
  "rnafusion.ConvertToCram_RNA.cram_exp": "*.cra[im]",
  "rnafusion.FastQC_RNA.html_exp": "*.html",
  "rnafusion.FastQC_RNA.zip_exp": "*.zip",
  "rnafusion.GatherQCreports_RNA.multiqc_command_line_config": "{ report_comment: This report was generated for the <code>RNA-Seq Gene Fusion and SNV</code> pipeline }",
  "rnafusion.MergeBamAlignment_RNA.bwa_commandline": null,
  "rnafusion.MergeBamAlignment_RNA.bwa_version": null,
  "rnafusion.MergeBamAlignment_RNA.program_group_name": null,
  "rnafusion.MergeBamAlignment_RNA.program_record": null,
  "rnafusion.StarAlignBam.alignIntronMax": "100000",
  "rnafusion.StarAlignBam.alignMatesGapMax": "100000",
  "rnafusion.StarAlignBam.alignSJDBoverhangMin": "10",
  "rnafusion.StarAlignBam.alignSJstitchMismatchNmax": null,
  "rnafusion.StarAlignBam.chimJunctionOverhangMin": "12",
  "rnafusion.StarAlignBam.chimSegmentMin": "12",
  "rnafusion.StarGenerateReferences.read_length": "150",
  "rnafusion.ValidateMarkDupSam_RNA.ignore": null,
  "rnafusion.ValidateMarkDupSam_RNA.max_output": null,
  "rnafusion.ValidateRecalibratedSam_RNA.max_output": null,
  "rnafusion.ValidateSam_RNA.ignore": null,
  "rnafusion.ValidateSam_RNA.max_output": null,
  "rnafusion.ValidateSplittedReadsSam_RNA.max_output": null,
  "rnafusion.biomaterial_id": "", #Adjust
  "rnafusion.blacklist": "/hpc/pmc_gen/annotation/genes/BlackList_fusion_genes_Homo_sapiens_GRCh38_GTF_v1.0.txt",
  "rnafusion.cache_path": "/hpc/pmc_kool/fvalzano/ubams/Some_sequencer/", #Adjust
  "rnafusion.cram_suffix": ".cram",
  "rnafusion.data_host": "hpct03",
  "rnafusion.dbSNP_vcf": "/hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
  "rnafusion.dbSNP_vcf_index": "/hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
  "rnafusion.derived_data_path": "",
  "rnafusion.experiment_ids": [
    "PMCRX000AJY"
  ],
  "rnafusion.experiment_type": "transcriptomics",
  "rnafusion.fc_report_reads_bam": "false",
  "rnafusion.fc_tmpspace_multiplier": "1.5",
  "rnafusion.features_to_quantify": [
    "gene",
    "exon"
  ],
  "rnafusion.genefusion_detection_rna_wf.min_ffpm": "0",
  "rnafusion.image_ensemblvep_digest": "docker.io/ensemblorg/ensembl-vep@sha256:62d189098d86cbee612ad761a2b7b88b664ebac36c5f77d10fc7c429e1d152ea",
  "rnafusion.image_fastqc_digest": "docker.io/biocontainers/fastqc@sha256:387748462c7fc280b7959ceda0f6251190d2e4b9ebc0585d24e7bcb58bdcf2bf",
  "rnafusion.image_gatk_digest": "docker.io/broadinstitute/gatk@sha256:f2602e0bbc0117c30d23d8d626eb8d0a21ca672bb71180b5cf25425603a0ae09",
  "rnafusion.image_htslib_digest": "docker.io/princessmaximacenter/htslib@sha256:893b7fbe342fba0a64277bd6051b2a008a8480d1167da6c217398a93b3871060",
  "rnafusion.image_multiqc_digest": "docker.io/princessmaximacenter/multiqc@sha256:f3e0ad418d11bef7d7db506ae4719787ce906d9e7c2569bbed84e042aea671a6",
  "rnafusion.image_picard_digest": "docker.io/princessmaximacenter/picard@sha256:d561453a4a25dfac7933ec2237259d1c5adc974287cc639418bd3c2af9881019",
  "rnafusion.image_python_digest": "docker.io/python@sha256:98fb5342195e69ffda54a7584ed202be71154c7ef64931da5bec5a41739c78d5",
  "rnafusion.image_rnaseqtools_digest": "docker.io/princessmaximacenter/rnaseq-tools@sha256:5ddf2ba2b1b95453dc3dd8f618595d67bb5fb2cfae44b6d3f7919038ed2c30a3",
  "rnafusion.image_snvtools_digest": "docker.io/princessmaximacenter/snv-tools@sha256:87619cc2ca94685696bdc57a7fab7f9af62697cf77c4d73df9509ba0d8019174",
  "rnafusion.image_starfusion_digest": "docker.io/sdevos/star-fusion@sha256:ee5c8ecb16e172db54d536b82917234865c1dde11722eb16327cc37e77469eed",
  "rnafusion.image_trecodetools_digest": "docker.io/princessmaximacenter/trecode-tools@sha256:dacf66d19ddf6fc81319d35c5aeab4620c588316c63e8a39647dfd395966dd9a",
  "rnafusion.input_files_md5": {
    "fileMd5": "", #Adjust with result of md5sum .ubam
    "fileUrl": "" #Adjust with ubam file
  },
  "rnafusion.known_indels_sites_VCFs": [
    "/hpc/pmc_gen/references/hg38bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf",
    "/hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.known_indels.vcf"
  ],
  "rnafusion.known_indels_sites_indices": [
    "/hpc/pmc_gen/references/hg38bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.idx",
    "/hpc/pmc_gen/references/hg38bundle/v0/Homo_sapiens_assembly38.known_indels.vcf.idx"
  ],
  "rnafusion.library_strategy": "RNA-Seq",
  "rnafusion.memory_huge": "48",
  "rnafusion.memory_moderate": "16",
  "rnafusion.memory_small": "8",
  "rnafusion.module_petalink_version": "",
  "rnafusion.module_samtools_version": "",
  "rnafusion.molgenis_token": "",
  "rnafusion.molgenis_url": "",
  "rnafusion.multiqc_config_file": "/hpc/pmc_gen/configurations/multiqc/rna_fusions_germline_snv_v1.0.1.yaml",
  "rnafusion.ref_ann": "/hpc/pmc_gen/annotation/genes/ref_annot_GRCh38_gencode_v31_CTAT_lib_Oct012019.gtf",
  "rnafusion.ref_dict": "/hpc/pmc_gen/references/RNA-Seq/references/ref_genome_GRCh38_gencode_v31_CTAT_lib_Oct012019.dict",
  "rnafusion.ref_fasta": "/hpc/pmc_gen/references/RNA-Seq/references/ref_genome_GRCh38_gencode_v31_CTAT_lib_Oct012019.fa",
  "rnafusion.ref_fasta_index": "/hpc/pmc_gen/references/RNA-Seq/references/ref_genome_GRCh38_gencode_v31_CTAT_lib_Oct012019.fa.fai",
  "rnafusion.ref_flat_ann": "/hpc/pmc_gen/annotation/genes/refFlat_GRCh38_gencode_v31_CTAT_lib_Oct012019_gtf.txt",
  "rnafusion.reference_version": "GRCh38",
  "rnafusion.ribosomal_interval_list": "/hpc/pmc_gen/annotation/genes/GRCh38_gencode_v31_CTAT_lib_Oct012019.ribosomalRNA.intervallist.txt",
  "rnafusion.run_ids": [
    "PMCRR000AJY"
  ],
  "rnafusion.scatter_count": "10",
  "rnafusion.snp_profile_target_bed_files": [
    "/hpc/pmc_gen/references/hg38bundle/v0/snp_profile_targets/PGX.bed.gz",
    "/hpc/pmc_gen/references/hg38bundle/v0/snp_profile_targets/WESWTS.bed.gz"
  ],
  "rnafusion.snp_profile_target_bed_index_files": [
    "/hpc/pmc_gen/references/hg38bundle/v0/snp_profile_targets/PGX.bed.gz.tbi",
    "/hpc/pmc_gen/references/hg38bundle/v0/snp_profile_targets/WESWTS.bed.gz.tbi"
  ],
  "rnafusion.starAlign_genomeDir": "GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/",
  "rnafusion.starFusion_genomeDir": "GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/",
  "rnafusion.store_primary_qa": false,
  "rnafusion.strip_url": "htt(p|ps)://[^/]*/",
  "rnafusion.study": [
    "PMCPR000AAF"
  ],
  "rnafusion.submission_date": "2021-09-20T15:43:49",
  "rnafusion.targetIntervalsBedFile": "/hpc/pmc_gen/references/hg38bundle/v0/genome_targets/hg38_diagnostic_somatic_2.1.bed.gz",
  "rnafusion.targetIntervalsBedIndexFile": "/hpc/pmc_gen/references/hg38bundle/v0/genome_targets/hg38_diagnostic_somatic_2.1.bed.gz.tbi",
  "rnafusion.target_host": "gerrit",
  "rnafusion.tmpspace_gb_huge": "80",
  "rnafusion.tmpspace_gb_moderate": "40",
  "rnafusion.tmpspace_gb_small": "5",
  "rnafusion.variant_filter": "'PASS','clustered_events'",
  "rnafusion.vep_assembly": "GRCh38",
  "rnafusion.vep_cache": "/hpc/pmc_gen/annotation/variants/cache/vep-104",
  "rnafusion.vep_fields": [
    "SYMBOL",
    "CCDS"
  ],
  "rnafusion.vep_species": "Homo sapiens",
  "rnafusion.do_virus_detection": true,
  "rnafusion.virus_genome_dir": "/hpc/pmc_gen/references/RNA-Seq/virus/virus_selection_v1.0.0.genomeDir.tar.gz",
  "rnafusion.virus_selection": "/hpc/pmc_gen/references/RNA-Seq/virus/virus_selection_v1.0.0.tsv",
  "rnafusion.wallclock_ages": "48:00:00",
  "rnafusion.wallclock_long": "24:00:00",
  "rnafusion.wallclock_moderate": "08:00:00",
  "rnafusion.wallclock_short": "04:00:00",
  "rnafusion.web_url": "",
  "rnafusion.wgs_calling_interval_list": "/hpc/pmc_gen/references/hg38bundle/v0/wgs_calling_regions.hg38.interval_list",
  "rnafusion.wgs_evaluation_interval_list": "/hpc/pmc_gen/references/hg38bundle/v0/qc/wgs_evaluation_regions.hg38.interval_list",
  "rnafusion.whitelist_geneids": "/hpc/pmc_gen/annotation/genes/WhiteList_fusion_genes_Homo_sapiens_GRCh38_BioMart_v1.7.txt",
  "rnafusion.workflow_instance_id": "PMCWI000AAN",
  "rnafusion.zipped_star_references": "/hpc/pmc_gen/references/RNA-Seq/starfusion/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play.tar.gz"
}'

        # Write the JSON content to the file
        echo "$json_content" > "$output_file"
    cd ..
done