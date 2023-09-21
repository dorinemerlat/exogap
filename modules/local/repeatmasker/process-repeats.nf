process PROCESS_REPEATS {
    tag "PROCESS_REPEATS_$genome_id"
    time '1d'
    publishDir "results/$genome_id/repeatmasker"

    conda (params.enable_conda ? 'py_fasta_validator==0.5--py39h7d875b9_0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/repeatmasker:4.1.5--pl5321hdfd78af_0':
        'quay.io/biocontainers/repeatmasker:4.1.5--pl5321hdfd78af_0' }"

    input:
    tuple val(masked_id), path(masked_path)
    tuple val(out_id), path(out_path)
    path catgz
    path repeats_lib

    output:
    tuple val(genome_id), path("${genome_id}.masked"),      emit: masked
    tuple val(genome_id), path("${genome_id}.gff"),         emit: gff
    tuple val(genome_id), path("${genome_id}-complex.gff"), emit: gff_complex
    tuple val(genome_id), path("${genome_id}.tbl"),         emit: tbl
    tuple val(genome_id), path("${genome_id}.out"),         emit: out
    tuple val(genome_id), path("${genome_id}.align"),       emit: align
    tuple val(genome_id), path("${genome_id}.cat"),         emit: cat

    script:
    """
    cat $catgz > ${genome_id}.cat
    cp $masked_path ${genome_id}.masked
    cp $out_path ${genome_id}.out

    ProcessRepeats -a -gff -lib $repeats_lib -xsmall $catgz

    # Isolation of complex repeats
    grep -v -e "Satellite" -e ")n" -e "-rich"  "${genome_id}.gff" > ${genome_id}-complex.gff.unreformated

    # Reformatting to make it work with Maker
    cat ${genome_id}-complex.gff.unreformated | \
        perl -ane '\$id; if(!/^\#/){{@F = split(/\t/, \$_); chomp \$F[-1];\$id++; \$F[-1] .= "\;ID=\$id"; \$_ =join("\t", @F)."\n"}} print \$_' \
        > ${genome_id}-complex.gff
    """
}
