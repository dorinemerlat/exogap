process GET_NEWICK {
    debug

    // conda (params.enable_conda ? 'bioconda requests==2.26.0' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/requests:2.26.0':
    //     'quay.io/biocontainers/requests:2.26.0' }"

    input:
    val(taxids)
    val(about_genomes)

    output:
    path('.newick')

    script:
    """
    taxids=(echo $taxids |sed "s/[][,]//g")
    about_genomes=(echo $about_genomes |sed "s/\], /\n/g" |sed -E "s/\[|\]//g" |sed "s/ //g" |tr '\n' ')
    get-newick.py -t "$taxids" -i "$about_genomes"
    """
}
