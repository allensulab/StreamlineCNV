#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.fastq = ''
params.species = 'Homo_sapiens'
params.assembly = 'GRCh38'
params.outdir = './results'
params.plotGeneDensity = false
params.geneList = ''  // GeneList file path for gene density plotting
params.recolor = false
params.sampleInfo = ''  // SampleInfo file path for recolor

fastq_ch = Channel.fromPath(params.fastq)
species_ch = Channel.fromPath(params.species)
assembly_ch = Channel.fromPath(params.assembly)

process TRIM {
    input:
    val inputfq

    output:
    path '*fq'

    script:
    sample_id = inputfq.name.split('\\.')[0]
    """
    fastx_trimmer -Q33 -i "$inputfq" -o "${sample_id}_trimmed.fq" -l 40
    """
}

process BWA {
    input:
    tuple(path(trimmedfq), path(species),path(assembly))

    output:
    path 'mapping_dir'

    script:
    sample_id = trimmedfq.name.split('_')[1]
    """
    bwa aln -t 1 /data/${species}.${assembly}.dna.primary_assembly.fa $trimmedfq >${sample_id}.sai
    bwa samse /data/${species}.${assembly}.dna.primary_assembly.fa ${sample_id}.sai $trimmedfq |grep -v MT|grep -v KI|grep -v GL|grep -v JH > ${sample_id}.sam
    samtools view -uSh ${sample_id}.sam | samtools sort - -T ${sample_id}.sort.tmp -o ${sample_id}.sort.bam
    samtools index ${sample_id}.sort.bam
    mkdir mapping_dir
    mv *bam* mapping_dir
    """
}

process HMMCOPY {
    publishDir params.outdir, mode: 'copy' 
    input:
    tuple(path(bam), path(assembly))

    output:
    path '*pdf', emit: figures_archive
    path '*tsv', emit: txt_results
    path '*wig', emit: wig
    script:
    sample_id = bam.name.split('_')[0]
    """
    python3 /opt/scripts/hmm_pipe.py --keeptmp $bam ${assembly} >stats.txt
    id=\$(ls *pdf | cut -f1 -d .)
    python3 /opt/scripts/segments_anno.py merge.segments.txt /data/${assembly} \${id}_merge.segments.tsv
    mv stats.txt \${id}_stat.tsv
    mv *tmp/input.wig ./\${id}.wig
    """
}

process GeneDensityPlot {
    when:
    params.plotGeneDensity

    input:
    tuple(path(geneList), path(assembly))

    output:
    path '*jpg', emit: GeneDensity_figures

    publishDir params.outdir, mode: 'copy', pattern: '*.jpg'

    script:
    """
    python3 /opt/scripts/RetrieveGeneLoc.py ${geneList} /data/${assembly}_chr Geneloc
    Rscript /opt/scripts/GeneDensityPlot.r ${assembly}
    """
}


process RECOLOR {
    when params.recolor

    input:
    tuple(path(sample_info_file), path(hmmcopy_wig_file))

    output:
    path '*pdf', emit: recolor_pdf
    publishDir params.outdir, mode: 'copy'

    script:
    // Extract the sample ID from the staged hmmcopy_wig_file name
    def id = hmmcopy_wig_file.getName().split("\\.")[0]

    // Subset sample_info_file to get the corresponding directory, assembly, and sex for this sample
    def sampleLine = file(sample_info_file).text.split("\\n").find { it.split("\\t")[1] == id }
    if (!sampleLine) {
        error "Sample ID ${id} not found in sample_info_file."
    }
    def fields = sampleLine.split("\\t")
    if (fields.size() < 3) {
        error "Sample info line for ${id} has less than 3 fields: ${sampleLine}"
    }
    def dir = fields[1]
    def assembly = fields[0]
    def sex = fields[2]

    """
    echo "Processing sample: ${id}"
    mkdir -p ${id}_recolor_output_dir
    cp ${hmmcopy_wig_file} ${id}_recolor_output_dir/input.wig

    if [ ! -f "${id}_recolor_output_dir/input.wig" ]; then
        echo "Error: ${id}_recolor_output_dir/input.wig not found."
        exit 1
    fi

    echo "Extracted assembly: ${assembly}"
    echo "Extracted sex: ${sex}"

    cp /opt/scripts/run_hmmcopy.${assembly}.${sex}.r ${id}_recolor_output_dir/
    cd ${id}_recolor_output_dir
    R CMD BATCH run_hmmcopy.${assembly}.${sex}.r

    # Move the output PDF to the mail folder
    mv single.fixed.recolor.pdf ../${id}_recolored.pdf
    """ 
}

workflow {
    trimresults_ch = TRIM(fastq_ch)
    bwatmp_ch = trimresults_ch.combine(species_ch).combine(assembly_ch)
    mappingresults_ch = BWA(bwatmp_ch)
    hmmtmp_ch = mappingresults_ch.flatten().combine(assembly_ch)
    hmmcopyresults_ch = HMMCOPY(hmmtmp_ch)

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

        hmm_wig_ch = hmmcopyresults_ch.wig
        log.info "### DEBUG: hmmcopyresults_ch.out.wig contents: ${hmm_wig_ch}"

        hmm_ch = sampleInfo_ch.combine(hmm_wig_ch)
        
        log.info "### DEBUG: Combined sampleInfo and wig channel contents: ${hmm_ch}"

        RECOLOR(hmm_ch)
    } else {
        log.info 'Skipping recolor as sampleInfo file is missing or recolor is false'
    }
}

