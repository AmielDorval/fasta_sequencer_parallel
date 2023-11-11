def fasta_to_list_str (uri) -> list[str]:
    seqList= []
    with open(uri, "r") as file:
        for line in file:
            if line[0] != ">":
                seqList.append(line.strip())
    return seqList