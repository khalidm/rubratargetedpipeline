
# stageDefaults contains the default options which are applied to each stage (command).
# This section is required for every Rubra pipeline.
# These can be overridden by options defined for individual stages, below.
# Stage options which Rubra will recognise are:
#  - distributed: a boolean determining whether the task should be submitted to a cluster
#      job scheduling system (True) or run on the system local to Rubra (False).
#  - walltime: for a distributed PBS job, gives the walltime requested from the job
#      queue system; the maximum allowed runtime. For local jobs has no effect.
#  - memInGB: for a distributed PBS job, gives the memory in Gigabytes requested from the
#      job queue system. For local jobs has no effect.
#  - queue: for a distributed PBS job, this is the name of the queue to submit the
#      job to. For local jobs has no effect.
#  - modules: the modules to be loaded before running the task. This is intended for
#      systems with environment modules installed. Rubra will call module load on each
#      required module before running the task. Note that defining modules for individual
#      stages will override (not add to) any modules listed here. This currently only
#      works for distributed jobs.
stageDefaults = {
    'distributed': True,
    'walltime': "08:00:00",
    'memInGB': 8,
    'queue': "batch",
    'jobscript': " --account VR0002",
    'modules': [
        "bwa-intel/0.7.5a",
        "samtools-intel/0.1.19",
        "picard/1.53",
        "python-gcc/2.7.5",
        "R-intel/2.15.3",
        "gatk/2.6-5"
    ]
}

pipeline = {
    'logDir': 'logfiles',
    'logFile': 'pipeline.log',
    'style': 'print',
    'procs': 20,
    'verbose': 1,
    'end': [ ],
    # 'end': [ 'earlyDepthOfCoverage', 'dedupedDepthOfCoverage', 'finalDepthOfCoverage',
    #         'fastqc',
    #         'igvcountMergedBams', 'countRunBam',
    #         'collateReadCounts',
    #         'vcfIndexSNPs', 'vcfIndexIndels',
    #         'getEnsemblAnnotations',
    #         'exonCoverage' ],
    'force': [],
    'restrict_samples': False
    # 'allowed_samples': ['060037']
    }

working_files = {
    'fastq_dirs': [
         # '/vlsci/VR0215/shared/khalid/targeted/testfastq'
         '/vlsci/VR0215/shared/khalid/targeted/fastq_files/kr'
    ],
    'fastq_symlink_dir': '/vlsci/VR0215/shared/khalid/targeted/outputkr/fastq_symlinks',
    'output_dir': '/vlsci/VR0215/shared/khalid/targeted/outputkr/',
    'varscanjar': '/vlsci/VR0002/kmahmood/Programs/varscan/VarScan.v2.3.6.jar',
    # 'snpeff': '/vlsci/VR0002/kmahmood/Programs/snpEff/'
    'snpeff': '/vlsci/VR0002/kmahmood/Programs/snpEff3.4/snpEff/'
}

ref_files = {
    'fasta_reference': '/vlsci/VR0215/shared/khalid/targeted/refs/genome.fa',
    'fasta_dict': '/vlsci/VR0215/shared/khalid/targeted/refs/genome.dict',
    # 'exon_bed': '/vlsci/VR0215/shared/khalid/targeted/target_beds/S0470822_Covered.bed',
    'exon_bed': '/vlsci/VR0215/shared/khalid/targeted/target_beds/S0470822_Covered_sorted_merged.bed',
    'exon_bed_extended': '/vlsci/VR0215/shared/khalid/targeted/target_beds/S0470822_Covered100.bed',
    'dbsnp': '/vlsci/VR0215/shared/khalid/targeted/refs/dbsnp_137.hg19.vcf',
    'indels_realign_goldstandard': '/vlsci/VR0215/shared/khalid/targeted/refs/Mills_and_1000G_gold_standard.indels.hg19.vcf',
    'indels_realign_1000G': '/vlsci/VR0215/shared/khalid/targeted/refs/1000G_phase1.indels.hg19.vcf'
    }

# stages should hold the details of each stage which can be called by runStageCheck.
# This section is required for every Rubra pipeline.
# Calling a stage in this way carries out checkpointing and, if desired, batch job
# submission.
# Each stage must contain a 'command' definition. See stageDefaults above for other
# allowable options.
stages = {
    "fastqc": {
        "command": "fastqc --quiet -o %outdir %seq",
        "walltime": "02:00:00",
        'modules': [ "fastqc/0.10.1" ]
    },
        'indexReferenceBWA': {
        'command': "bwa index %ref -a bwtsw",
        'walltime': "05:00:00"
    },
    'indexReferenceSAM': {
        'command': "samtools faidx %ref"
    },
    'alignBWA': {
        'command': "bwa aln -t 8 %encodingflag %ref %seq > %out",
        'walltime': "04:00:00",
        'procs': 8,
        # 'queue': 'smp',
        'memInGB': 23
    },
    'alignToSamSE': {
        'command': "bwa samse %ref %meta %align %seq > %out"
    },
    'alignToSamPE': {
        'command': "bwa sampe %ref %meta %align1 %align2 %seq1 %seq2 > %out",
        'memInGB': 32
    },
    'samToSortedBam': {
        'command': "./SortSam 23 VALIDATION_STRINGENCY=LENIENT INPUT=%seq OUTPUT=%out SORT_ORDER=coordinate",
        'memInGB': 32,
        'walltime': "05:00:00"
    },
    'mergeBams': {
        'command': "./PicardMerge 24 %baminputs USE_THREADING=true VALIDATION_STRINGENCY=LENIENT AS=true OUTPUT=%out",
        'procs': 2,
        'memInGB': 32,
        'walltime': "10:00:00"
    },
    'indexBam': {
        'command': "samtools index %bam"
    },
    'flagstat': {
        'command': "samtools flagstat %bam > %out",
        'walltime': "00:10:00"
    },
    'igvcount': {
        'command': "igvtools count %bam %out hg19",
        'modules': [ "igv/2.3.15" ]
    },
    'indexVCF': {
        'command': "./vcftools_prepare.sh %vcf",
        'modules': [ "tabix/0.2.5" ]
    },
    'dedup': {
        'command': "./MarkDuplicates 24 INPUT=%bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE=%log OUTPUT=%out",
        'memInGB': 32,
        'walltime': '05:00:00'
    },
    'realignIntervals': {
        # Hard-coded to take 2 known indels files right now
        'command': "./GenomeAnalysisTK 23 -T RealignerTargetCreator -R %ref -I %bam --known %indels_goldstandard --known %indels_1000G -L %bed -log %log -o %out",
        'memInGB': 32,
        'walltime': "08:00:00"
    },
    'realign': {
        'command': "./GenomeAnalysisTK 22 -T IndelRealigner -R %ref -I %bam -known %indels_goldstandard -known %indels_1000G -targetIntervals %intervals -log %log -o %out",
        'memInGB': 32,
        'walltime': "08:00:00"
    },
    # Original: moved above before realignemnt as recommended by GATK
    # 'dedup': {
    #     'command': "./MarkDuplicates 24 INPUT=%bam REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT AS=true METRICS_FILE=%log OUTPUT=%out",
    #     'memInGB': 32,
    #     'walltime': '05:00:00'
    # },
#GATK1    'baseQualRecalCount': {
#GATK1        'command': "./GenomeAnalysisTK 12 -T CountCovariates -I %bam -R %ref --knownSites %dbsnp -nt 8 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -log %log -recalFile %out",
#GATK1        'queue': 'smp',
#GATK1        'memInGB': 23,
#GATK1        'walltime': "3:00:00:00"
#GATK1    },
#GATK1    'baseQualRecalTabulate': {
#GATK1        'command': "./GenomeAnalysisTK 4 -T TableRecalibration -I %bam -R %ref -recalFile %csvfile -l INFO -log %log -o %out",
#GATK1        'walltime': "3:00:00:00"
#GATK1    },
    'leftalignindels': {
        'command': "./GenomeAnalysisTK 24 -allowPotentiallyMisencodedQuals -T LeftAlignIndels -I %input -R %ref -o %output",
        'memInGB': 32,
        'walltime': "05:00:00"
    },
    'baseQualRecal': {
        'command': "./GenomeAnalysisTK 24 -T BaseRecalibrator -I %bam -R %ref --knownSites %dbsnp --knownSites %indels_goldstandard -log %log -o %out",
        'memInGB': 32,
        'walltime': "05:00:00"
    },
    'baseQualRecalPrintReads': {
        'command': "./GenomeAnalysisTK 32 -T PrintReads -I %bam -R %ref -BQSR %csvfile -log %log  -o %out",
        'memInGB': 48,
        'walltime': "05:00:00"
    },
    'callSNPs': {
        'command': "./GenomeAnalysisTK 24 -T UnifiedGenotyper -nt 8 -R %ref -I %bam -L %bed --dbsnp %dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -l INFO -A AlleleBalance -A Coverage -A FisherStrand -glm SNP -log %log -o %out",
        # 'queue': 'smp',
        'procs': 8,
        'memInGB': 32,
        'walltime': "10:00:00"
    },
    'callHAP': {
        'command': "./GenomeAnalysisTK 24 -T HaplotypeCaller -R %ref -I %bam -L %bed --dbsnp %dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -l INFO -A AlleleBalance -A Coverage -A FisherStrand -log %log -o %out",
        # 'queue': 'smp',
        'procs': 2,
        'memInGB': 32,
        'walltime': "12:00:00"
    },
    # 'callHAPMerged': {
    #     'command': "./GenomeAnalysisTK 24 -T HaplotypeCaller -R %ref %bam -L %bed --dbsnp %dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -l INFO -A AlleleBalance -A Coverage -A FisherStrand -log %log -o %out",
    #     'procs': 2,
    #     'memInGB': 32,
    #     'walltime': "240:00:00"
    # },
    'callIndels': {
        'command': "./GenomeAnalysisTK 24 -T UnifiedGenotyper -nt 8 -R %ref -I %bam -L %bed --dbsnp %dbsnp -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200 -l INFO -A AlleleBalance -A Coverage -A FisherStrand -glm INDEL -log %log -o %out",
        'procs': 8,
        'memInGB': 32,
        'walltime': "12:00:00"
    },
    'callVARSCAN': {
        'command': "./RunVarscan %ref %varscanjar %samplelist %out %bam",
        'procs': 1,
        'memInGB': 32,
        'modules': [ "samtools-intel/0.1.19", "java/1.7.0_25" ],
        'walltime': "8:00:00"
    },
    # 'callVARSCANMerged': {
    #     'command': "./RunVarscan %ref %varscanjar %samplelist %out %bam",
    #     # 'command': "./RunVarscan %ref %bam %varscanjar %out",
    #     'procs': 1,
    #     'memInGB': 32,
    #     'modules': [ "samtools-intel/0.1.19", "java/1.7.0_25" ],
    #     'walltime': "48:00:00"
    # },
    'filterSNPs': {
        # Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
        'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
        'memInGB': 32,
    },
    'filterVarscan': {
        # Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
        #'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
        'command': "./RunVarscanFilter %ref %varscanjar %vcf %out",
        'procs': 1,
        'memInGB': 32,
        'modules': [ "samtools-intel/0.1.19", "java/1.7.0_25" ]
    },
    # 'filterVarscan2': {
    #     # Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
    #     # 'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
    #     'command': "./RunVarscan %ref %varscanjar %vcf %out",
    #     'memInGB': 32,
    # },
    'filterHapVcfs': {
        'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
        'memInGB': 32,
    },
    # 'filterMergedHapVcfs': {
    #     'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
    #     'memInGB': 32,
    # },
    'filterIndels': {
        # Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
        # If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
        'command': "./GenomeAnalysisTK 4 -T VariantFiltration -R %ref --variant %vcf --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0' --filterName 'GATK_MINIMAL_FILTER' -log %log -o %out",
        'memInGB': 32,
    },
    # 'annotateEnsembl': {
    #     # This command as written assumes that VEP and its cache have been
    #     # downloaded in respective locations
    #     # ./variant_effect_predictor_2.5
    #     # ./variant_effect_predictor_2.5/vep_cache
    #     'command': "perl variant_effect_predictor_2.5/variant_effect_predictor.pl --cache --dir variant_effect_predictor_2.5/vep_cache -i %vcf --vcf -o %out -species human --canonical --gene --protein --sift=b --polyphen=b > %log",
    #     'modules': [ "perl/5.18.0", "ensembl/67" ]
    # },
     'annotateSNPEff': {
        # This command as written assumes that snpEFF database hg19 is downloaded and the config file is amended accordingly
        'command': "./SnpEff 24 %snpeff eff -c %config -v hg19 %vcf > %output",
        'memInGB': 32,
        'walltime': "12:00:00"
    },
    'depthOfCoverage': {
        'command': "./GenomeAnalysisTK 24 -T DepthOfCoverage -R %ref -I %bam -L %bed -omitBaseOutput -ct 1 -ct 10 -ct 20 -ct 30 -o %out",
        'memInGB': 32,
    },
    'exonCoverage': {
        'command': "coverageBed -abam %bam -b %exon_bed > %out",
        'modules': [ "bedtools-intel/2.17.0" ]
    },
    'intersectBam': {
        'command': "intersectBed -abam %bam -b %bed > %out",
        'modules': [ "bedtools-intel/2.17.0" ]
    },
    'collateReadcounts': {
        'command': 'python count_flagstat_exome.py %dir %outdir',
        'walltime': "00:10:00"
    }
}
