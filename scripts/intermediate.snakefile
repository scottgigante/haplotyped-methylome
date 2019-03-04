rule download_nanopolish_b6xcast:
    output:
        temp("../nanopore/b6xcast.minion.methylation.tsv.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERZ784/ERZ784839/albacore_1.2.2.b6xcast.methylation.tsv.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_nanopolish_castxb6:
    output:
        temp("../nanopore/castxb6.promethion.methylation.tsv.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERZ784/ERZ784841/albacore_2.2.7.castxb6.promethion.methylation.tsv.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_nanopolish_cast:
    output:
        temp("../nanopore/cast.minion.methylation.tsv.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERZ784/ERZ784840/albacore_1.2.2.cast.methylation.tsv.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule download_nanopolish_b6:
    output:
        temp("../nanopore/b6.minion.methylation.tsv.gz")
    params:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/ERZ784/ERZ784838/albacore_1.2.2.b6.methylation.tsv.gz",
    shell:
        "wget -q -O {output} {params.url}"

rule intermediate_md5:
    input:
        tsv = "../nanopore/{archive}.methylation.tsv.gz",
        md5 = ancient("../md5/nanopore/{archive}.methylation.tsv.gz.md5")
    output:
        temp("../nanopore/{archive}.methylation.tsv.gz.md5_ok")
    shell:
        "md5sum -c {input.md5} && touch {output}"

rule intermediate_gunzip:
    input:
        tsv = "../nanopore/{archive}.methylation.tsv.gz",
        md5_ok = "../nanopore/{archive}.methylation.tsv.gz.md5_ok"
    output:
        signpost = "../nanopore/{archive}.intermediate_download"
    shell:
        "gunzip {input.tsv} && touch {output.signpost}"
