process GALBA {
    tag "${id}"
    cpus 10
    time '10d'
    label 'braker'
    scratch false

    input:
    tuple val(id), val(meta), path(genome), path(proteins), val(specie_name)

    output:
    tuple val(id), val(meta), path("galba_${id}"), emit: galba_dir

    script:
    def name    = specie_name.replaceAll(/\s+/, '_').replace('.','')
    def toks    = name.tokenize('_')
    def genus   = toks.size() > 0 ? toks[0] : name
    def species = toks.size() > 1 ? toks[1..-1].join('_') : ''

    """
    cp -r /opt/Augustus/config/ .
    AUGUSTUS_CONFIG_PATH="\$PWD/config"

    galba.pl --species=${genus}_${species} --genome=${genome} --prot_seq=${proteins} --gff3
    """

    stub:
    """
    touch galba_${id}.gff3
    """
}
