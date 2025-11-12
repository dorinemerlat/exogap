process CLEAN_REPEAT_FAMILIES {
    tag "${id}"
    label 'hite'
    scratch 'false'

    input:
    tuple val(id), val(meta), path(families), path(classification), val(source)

    output:
    tuple val(id), val(meta), path("${families.baseName}.cleaned.fa"), emit: main
    tuple val(id), val(meta), path("${families.baseName}.cleaned_for_mchelper.fa"), emit: for_mchelper

    script:
    """
    clean_repeat_families.py -f ${families} -c ${classification} -s ${source} -m ${meta.mnemonic.exogap}
    """

    stub:
    """
    touch ${families.baseName}.cleaned.fa ${families.baseName}.cleaned_for_mchelper.fa
    """
}
