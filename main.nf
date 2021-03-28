nextflow.enable.dsl = 2


process assembly {
    container 'nanozoo/shovill:1.1.0--1dafaa5'
    publishDir "${params.results}", mode: 'copy', overwrite: true
    cpus = 8
    memory '14 GB'

    input:
        tuple(val(name), path(genomes))

    output:
        tuple(val(name), path("*.fasta"))

    """
    shovill --R1 ${genomes[0]} --R2 ${genomes[1]} --gsize 5M --assembler megahit --trim --outdir assembly --minlen ${params.minlen} --cpus ${task.cpus} --force
    mv assembly/contigs.fa ${name}.fasta
    """
}


process annotate {
    container 'nanozoo/prokka:1.14.6--773a90d'
    publishDir "${params.results}", mode: 'copy', overwrite: true
    cpus = 8

    input:
        tuple(val(name), path(genomes))

    output:
        tuple(val(name), path("${name}.gff"))
        tuple(val(name), path("${name}.faa"), emit: 'proteins')

    """
    prokka --mincontiglen ${params.minlen} --cpus ${task.cpus} --outdir anno --prefix ${name} ${genomes}
    cp anno/${name}.gff ${name}.gff
    cp anno/${name}.faa ${name}.faa
    """
}


process checkm {
    /*
    Adapted from Muffin wf:

    https://github.com/RVanDamme/MUFFIN/blob/master/modules/checkm.nf
    */
    container 'nanozoo/checkm:1.1.3--c79a047'
    //maxForks 1
    publishDir "${params.results}/checkm/${name}/", mode: 'copy', pattern: "summary.txt"
    publishDir "${params.results}/checkm/${name}/", mode: 'copy', pattern: "taxonomy.txt"
    publishDir "${params.results}/checkm/${name}/", mode: 'copy', pattern: "*_checkm"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    cpus = 16
    memory '20 GB'

    input:
        tuple val(name), path(assembly)
    
    output:
        tuple val(name), path("summary.txt")
        tuple path("${name}_checkm"), path("taxonomy.txt")
    
    """
    mkdir tmp
    mkdir input
    cp ${assembly} input/assembly.fa
    checkm lineage_wf --tmpdir tmp --pplacer_threads 4 -t ${task.cpus} --reduced_tree -x fa input ${name}_checkm > summary.txt
    checkm tree_qa ${name}_checkm > taxonomy.txt
    """
}


process concern {
    container 'nanozoo/mmseqs2:11.e1a1c--55acb62'
    publishDir "${params.results}", mode: 'copy', overwrite: true
    // TODO: Create db beforehand
    
    input:
        tuple(val(name), path(proteins), path(db))

    output:
        tuple(val(name), path('aln.m8'))

    """
    mmseqs easy-search --max-accept 1 --min-seq-id 0.8 -c 0.5 ${proteins} ${db} aln.m8 tmp
    """
}


workflow {
    if (params.sra) {
        reads = channel.fromSRA(params.sra)
    }
    else {
        reads = channel.fromPath(params.genomes, checkIfExists: true)
                       .splitCsv(header: false)
                       .map{ row -> tuple(row[0], [row[1], row[2]]) }
    }

    assembly(reads)
    annotate(assembly.out)

    if (params.qc) {
        checkm(assembly.out)    
    }
    
    db = channel.fromPath(params.db, checkIfExists: true)
    concern(annotate.out.proteins.combine(db))
    // TODO: review(concern.out)

}
