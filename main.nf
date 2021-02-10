nextflow.enable.dsl = 2


process assembly {
    container 'nanozoo/spades:3.14.1--794e6a2'
    publishDir "${params.results}", mode: 'copy', overwrite: true
    cpus = 8

    input:
        tuple(val(name), path(genomes))

    output:
        tuple(val(name), path("*.fasta"))

    """
    shovill --R1 ${genomes[0]} --R2 ${genomes[1]} --gsize 5M --assembler spades --trim --outdir assembly --minlen ${params.minlen} --cpus ${task.cpus} --force
    mv assembly/contigs.fa ${name}.fasta
    """
}


process annotate {
    container 'nanozoo/prokka:1.13.4--d6a71cb'
    publishDir "${params.results}", mode: 'copy', overwrite: true
    cpus = 8

    input:
        tuple(val(name), path(genomes))

    output:
        tuple(val(name), path("${name}.gff3"))

    """
    prokka --mincontiglen ${params.minlen} --cpus ${task.cpus} --outdir anno --prefix ${name} 
    mv anno/${name}.gff3 ${name}.gff3
    """

}

process checkm {
    /*
    Adapted from Muffin wf:

    https://github.com/RVanDamme/MUFFIN/blob/master/modules/checkm.nf
    */
    container 'nanozoo/checkm:1.1.3--c79a047'
    maxForks 1
    publishDir "${params.results}/checkm/${name}/", mode: 'copy', pattern: "summary.txt"
    publishDir "${params.results}/checkm/${name}/", mode: 'copy', pattern: "taxonomy.txt"
    publishDir "${params.results}/checkm/${name}/", mode: 'copy', pattern: "*_checkm"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    
    input:
        tuple val(name), path(assembly)
    
    output:
        tuple val(name), path("summary.txt")
        tuple path("${name}_checkm"), path("taxonomy.txt")
    
    """
    mkdir tmp
    mkdir input
    mv ${assembly} input/assembly.fa
    checkm lineage_wf --tmpdir tmp --pplacer_threads 4 -t ${task.cpus} --reduced_tree -x fa input ${name}_checkm > summary.txt
    checkm tree_qa ${name}_checkm > taxonomy.txt
    """
}


workflow {
    if (params.sra) {
        reads = channel.fromSRA("${params.sra}")
    }
    else {
        reads = channel.fromPath(params.genomes, checkIfExists: true)
                       .splitCsv(header: false)
                       .map{ row -> tuple(row[0], row[1], row[2]) }
    }

    assembly(reads)
    annotate(assembly.out)
    checkm(assembly.out)
}
