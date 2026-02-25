#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.fastq = ''
params.species = 'Homo_sapiens'
params.assembly = 'GRCh38'
params.outdir = './results'
params.plotGeneDensity = false
params.clustering = false
params.clusteringLabel = ''  // Clustering label file path for sample clustering
params.geneList = ''  // GeneList file path for gene density plotting
params.recolor = false
params.sampleInfo = ''  // SampleInfo file path for recolor
params.clustering_memory = '16.GB'
params.no_sampleLabel = false
params.clustering_cpus = 4
params.clustering_time = '10.h'
params.tissue_color_palette = 'Set3'
params.dropChr = '' // space-separated list of chromosomes to drop, e.g. 'X Y MT'

fastq_ch = Channel.fromPath(params.fastq, type: 'file')
species_ch = Channel.fromPath(params.species)
assembly_ch = Channel.fromPath(params.assembly)
clusteringLabel_ch = Channel.fromPath(params.clusteringLabel, type: 'file')

process FASTQC {
    tag { inputfq.baseName }

    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
    path inputfq

    output:
    path "fastqc_*_logs"

    script:
    sample_id = inputfq.baseName.replaceAll(/\.gz$|\.fq$|\.fastq$/, '')
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q "$inputfq"
    """
}


process CHECK_FASTQ {
    tag "$inputfq"

    publishDir "${params.outdir}/QC", mode: 'copy', pattern: '*warning.txt'

    input:
    path inputfq

    output:
    tuple path("*.checked.fq"), val(true), emit: valid_fastq, optional: true
    tuple file("*_warning.txt"), val(true), optional: true
    tuple file("*int"), val(true), optional: true

    script:
    sample_id = inputfq.baseName.replaceAll(/\.gz$|\.fq$|\.fastq$/, '')
    """
    fqfile="${inputfq}"

    if [[ "${inputfq}" == *.gz ]]; then
        gunzip -c "${inputfq}" > "${sample_id}.int"
        fqfile="${sample_id}.int"
    fi

    nseq=\$(cat "\$fqfile" | wc -l)
    nreads=\$((nseq / 4))

    if [ "\$nreads" -lt 50000 ]; then
        echo "WARNING: ${inputfq} contains only \$nreads reads. Skipping. It is insufficient for reliable CNV estimation." > "${sample_id}_warning.txt"
    else
        cp "\$fqfile" "${sample_id}.checked.fq"
        touch "${sample_id}_passed_check.txt"
    fi
    """
}

process TRIM {
    input:
    tuple path(inputfq), val(valid)

    output:
    path '*fq'

    script:
    sample_id = inputfq.name.split('\\.')[0]
    """
    fastx_trimmer -Q33 -i "$inputfq" -o "${sample_id}_trimmed.fq" -l 40
    """
}

process BWA {
    publishDir "${params.outdir}/QC", mode: 'copy', pattern: '*.{report}' 
    input:
    tuple(path(trimmedfq), path(species),path(assembly))

    output:
    path 'mapping_dir', emit: mappingdir
    path '*report', emit: QC_report

    script:
    sample_id = trimmedfq.name.split('_')[1]
    """
    bwa aln -t 1 /data/${species}.${assembly}.dna.primary_assembly.fa $trimmedfq >${sample_id}.sai
    bwa samse /data/${species}.${assembly}.dna.primary_assembly.fa ${sample_id}.sai $trimmedfq |grep -v MT|grep -v KI|grep -v GL|grep -v JH > ${sample_id}.sam
    samtools view -uSh ${sample_id}.sam | samtools sort - -T ${sample_id}.sort.tmp -o ${sample_id}.sort.bam
    samtools index ${sample_id}.sort.bam
    denom=\$(tail -n1 /data/${assembly}_chrsize | cut -f2)
    echo "Sequencing coverage after trimming:" >>${sample_id}_sequencing_coverage.report
    samtools view -c -F 4 ${sample_id}.sort.bam | awk -v denom=\$denom '{printf "%.6f\\n", \$1 * 40 / denom}' >> ${sample_id}_sequencing_coverage.report
    mkdir mapping_dir
    mv *bam* mapping_dir
    """
}

process HMMCOPY {
    publishDir "${params.outdir}/HMMCOPY_results", mode: 'copy', pattern: '*.{pdf,txt,tsv,wig}' 
    input:
    tuple(path(bam), path(assembly))

    output:
    path '*pdf', emit: figures_archive
    path '*txt', emit: txt_results
    path '*tsv', emit: tsv_results
    path '*wig', emit: wig
    script:
    sample_id = bam.name.split('_')[0]
    """
    /usr/bin/python3 /opt/scripts/hmm_pipe.py --keeptmp $bam ${assembly} >stats.txt
    id=\$(ls *pdf | cut -f1 -d .)
    /usr/bin/python3 /opt/scripts/segments_anno.py merge.segments.txt /data/${assembly} \${id}_merge.segments.tsv
    /usr/bin/python3 /opt/scripts/aneuploidySum.py \${id}_merge.segments.tsv \${id}.int \${id}
    /usr/bin/python3 /opt/scripts/cnvScore.py /data/${assembly}_chrsize \${id}.int \${id}_cnvScore.tsv
    tail -n2 stats.txt > \${id}_stat.txt
    rm stats.txt
    mv *tmp/input.wig ./\${id}.wig
    """
}

process CNVSCORE {
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.{txt}'

    input:
    path segment_files


    output:
    path "Merge_seg_CNVscore.txt", emit: merged_segments

    script:
    // Prepare bash-safe list of file names for debug message

    """
    head -n1 ${baseDir}/${params.outdir}/HMMCOPY_results/*_cnvScore.tsv |sort -u |tail -n1> header
    echo "Reading from: ${baseDir}/${params.outdir}/HMMCOPY_results/"
    cat ${baseDir}/${params.outdir}/HMMCOPY_results/*_cnvScore.tsv | grep -v sampleID > int
    cat header int >Merge_seg_CNVscore.txt

    """
}

process CLUSTERING {
    cpus   params.clustering_cpus
    memory params.clustering_memory

    tag "${label_file}"
    publishDir "${params.outdir}/CLUSTERING_results", mode: 'copy', pattern: '*.{pdf}'

    input:
    path label_file
    path segment_files
    path assembly

    output:
    path "seg.txt", emit: merged_segments
    path "header", emit: header_file, optional: true
    path "seg_anno", emit: annotated_segments, optional: true
    path "clustering.pdf", emit: clusteringpdf, optional: true

    script:
    // Optional flag for Python
    def noSampleLabelOpt = params.no_sampleLabel ? '--no-sampleLabel True' : ''
    def dropChrOpt = params.dropChr ? "--drop-chr ${params.dropChr}" : ''
    def tissuePaletteOpt = params.tissue_color_palette ? "--palette ${params.tissue_color_palette}" : ''

    """
    echo "=== DEBUG INFO ==="
    echo "Current working directory: \$(pwd)"
    echo "Available files:"
    ls -la
    
    head -n 1 ${label_file} > header
    echo "Reading from: ${baseDir}/${params.outdir}/HMMCOPY_results/"
    echo "Label file: ${label_file}"
    cat ${baseDir}/${params.outdir}/HMMCOPY_results/*_merge.segments.tsv | grep -v SampleName | sed -e 's/\\.sort.tmp\\/segments.txt//g' | sed 's/"//g' | cut -f1,3-6 > seg.txt
    /usr/bin/python3 /opt/scripts/overlap.py seg.txt ${label_file} seg_anno
    /usr/bin/python3 /opt/scripts/clustering.py -c /data/${assembly}_chrsize -i seg_anno -o clustering.pdf --bin-size 500000 --neutral 3 --anchor 0 --coords one_based_inclusive ${noSampleLabelOpt} ${dropChrOpt} ${tissuePaletteOpt}
    echo "Final files created:"
    ls -la
    
    """
}



process GeneDensityPlot {
    when:
    params.plotGeneDensity

    input:
    tuple(path(geneList), path(assembly))

    output:
    file "*.jpg"

    publishDir "${params.outdir}/GeneDensityPlot", mode: 'copy', pattern: "*.jpg"

    script:
    """
    /usr/bin/python3 /opt/scripts/RetrieveGeneLoc.py ${geneList} /data/${assembly}_chr Geneloc
    Rscript /opt/scripts/GeneDensityPlot.r ${assembly}
    """
}


process RECOLOR_SAMPLES {
    publishDir "${params.outdir}/recolor", mode: 'copy', pattern: 'warning.txt'

    input:
    path sample_info_file
    path hmmcopy_wig_file

    output:
    path "clean_sampleInfo.txt", emit: cleanSampleInfo
    path "warning.txt", emit: warning_file, optional: true

    script:

    """
    ls ${baseDir}/${params.outdir}/HMMCOPY_results/*wig >wig
    python /net/ostrom/data/bcc/projects/allen-Su/2025/SL-CNV/bin/clean_sampleInfo.py ${sample_info_file} wig clean_sampleInfo.txt warning.txt
    if [ -f warning.txt ] && [ -s warning.txt ]; then
        echo "Warning file contains content - will be published"
    else
        rm -f warning.txt
        echo "Warning file is empty or doesn't exist - removed"
    fi
    
    """
}

process RECOLOR {
    publishDir "${params.outdir}/recolor", mode: 'copy', pattern: '*.{pdf,txt}'
    when: params.recolor

    input:
    tuple path(sample_info_file), path(hmmcopy_wig_file)

    output:
    path '*pdf', emit: recolor_pdf, optional: true

    script:
    def id = hmmcopy_wig_file.getName().split("\\.")[0]
    """
    echo "Processing sample: ${id}"
    echo "Sample info file: ${sample_info_file}"
    echo "Wig file: ${hmmcopy_wig_file}"

    ls -la ${sample_info_file}
    ls -la ${hmmcopy_wig_file}

    sample_line=\$(grep -P "^[^\\t]*\\t${id}\\t" ${sample_info_file} || echo "")
    
    assembly=\$(echo "\$sample_line" | cut -f1)
    sample_id=\$(echo "\$sample_line" | cut -f2)
    sex=\$(echo "\$sample_line" | cut -f3)
    
    echo "Found sample info - Assembly: \$assembly, Sample: \$sample_id, Sex: \$sex"

    if [ -z "\$assembly" ] || [ -z "\$sample_id" ] || [ -z "\$sex" ]; then
        echo "ERROR: Incomplete sample info for ${id}: \$sample_line" >> ${id}_warning.txt
    fi
    
    mkdir -p ${id}_recolor_output_dir
    cp ${hmmcopy_wig_file} ${id}_recolor_output_dir/input.wig
    r_script="/opt/scripts/run_hmmcopy.\${assembly}.\${sex}.r"

    cp "\$r_script" ${id}_recolor_output_dir/
    cd ${id}_recolor_output_dir
    
    if R CMD BATCH run_hmmcopy.\${assembly}.\${sex}.r; then
        if [ -f "single.fixed.recolor.pdf" ]; then
            mv single.fixed.recolor.pdf ../${id}_recolored.pdf
            echo "Successfully processed sample: ${id}"
        else
            echo "ERROR: R script did not generate PDF" >>${id}_warning.txt
        fi
    else
        echo "ERROR: R script execution failed" >>${id}_warning.txt
    fi
    """
}


workflow {
    fastqc_ch = FASTQC(fastq_ch)
    fastq_checked_ch = fastq_ch | CHECK_FASTQ
    trimresults_ch = TRIM(fastq_checked_ch.valid_fastq)
    bwatmp_ch = trimresults_ch.combine(species_ch).combine(assembly_ch)
    mappingresults_ch = BWA(bwatmp_ch)
    hmmtmp_ch = mappingresults_ch.mappingdir.flatten().combine(assembly_ch)
    hmmcopyresults_ch = HMMCOPY(hmmtmp_ch)
    cnvscore_input_ch = hmmcopyresults_ch.tsv_results.collect()
    cnvscore_ch = CNVSCORE(cnvscore_input_ch)

    if (params.clustering) {
        clustering_Label_ch = Channel.fromPath(params.clusteringLabel)
        segment_files_ch = hmmcopyresults_ch.tsv_results.collect()
        CLUSTERING(clustering_Label_ch, segment_files_ch,assembly_ch)
    }


    if (params.plotGeneDensity) {
        Channel
            .fromPath(params.geneList)
            .ifEmpty { error 'GeneList file is required for gene density plotting!' }
            .set { geneList_ch }
        geneListtmp_ch = geneList_ch.combine(assembly_ch)
        GeneDensityPlot(geneListtmp_ch)
    }

    if (params.recolor && params.sampleInfo) {
        sampleInfo_ch = Channel
            .fromPath(params.sampleInfo, checkIfExists: true)
            .ifEmpty { error 'sampleInfo file is required for recolor!' }
        
        log.info "### DEBUG: sampleInfo_ch contents: ${sampleInfo_ch}"
        wig_files_ch = hmmcopyresults_ch.wig.collect()
        hmm_wig_ch = hmmcopyresults_ch.wig
        log.info "### DEBUG: hmmcopyresults_ch.out.wig contents: ${hmm_wig_ch}"

        cleansample_ch=RECOLOR_SAMPLES(sampleInfo_ch,wig_files_ch)

        hmm_ch = cleansample_ch.cleanSampleInfo.combine(hmm_wig_ch)
        
        log.info "### DEBUG: Combined sampleInfo and wig channel contents: ${hmm_ch}"

        RECOLOR(hmm_ch)
    } else {
        log.info 'Skipping recolor as sampleInfo file is missing or recolor is false'
    }
}

