# Count CpG in a fasta file
# Xuanken Tay, 2018

import sys
import error
import re


def count_cpg(file, path=True):
    substring = "CG"
    print("chr\tstart\tend")
    f = open(file) if path else file
    line = f.readline()
    while line:
        if line.startswith('>'):
            chrm = line[1:].rstrip()
        else:
            index = 0
            linelen = len(line)
            # Find first non-N base
            # index = re.search(r'[^N]', line).start()
            # print("start\t{}\tend\t{}".format(index+1, linelen))
            while index < linelen:
                index = line.find(substring, index)
                if index == -1:
                    break
                print("{}\t{}\t{}".format(
                    chrm, index + 1, index + len(substring)))
                index += len(substring)
        line = f.readline()
    if path:
        f.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        if sys.stdin:
            count_cpg(sys.stdin, path=False)
        else:
            error.print_error()
    else:
        count_cpg(sys.argv[1])
