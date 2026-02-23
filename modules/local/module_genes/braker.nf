process BRAKER {
    tag "${id}"
    cpus 10
    time '10d'
    label 'braker'
    scratch false

    input:
    tuple val(id), val(meta), path(genome), path(bam), path(proteins), val(specie_name)

    output:
    tuple val(id), val(meta), path("braker_${id}"), emit: braker_dir

    script:
    def name    = specie_name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''

    """
    cp -r /opt/Augustus/config/ .
    AUGUSTUS_CONFIG_PATH="\$PWD/config"

    # if bam file is not empty
    if [[ -s ${bam} ]]; then
        bam_arg="--bam=${bam}"
    else
        bam_arg=""
    fi

    braker.pl --genome=${genome} --prot_seq=${proteins} \$bam_arg --gff3 --species=${genus}_${species} --threads ${task.cpus} \
        --AUGUSTUS_CONFIG_PATH=\$AUGUSTUS_CONFIG_PATH \
                                --AUGUSTUS_BIN_PATH="/opt/Augustus/bin/" \
                                --AUGUSTUS_SCRIPTS_PATH="/opt/Augustus/scripts/" \
                                --DIAMOND_PATH="/opt/ETP/tools/" \
                                --PROTHINT_PATH="/opt/ETP/bin/gmes/ProtHint/bin" \
                            --GENEMARK_PATH=/opt/ETP/bin/   --workingdir=./braker_${id}
    """

    stub:
    """
    touch braker.gff3
    """
}
