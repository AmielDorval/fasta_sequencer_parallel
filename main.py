import fastaReader
import timeit
import multiprocessing
from functools import partial


def best_overlap_with_min(a: str, b: str, min_overlap: int = 0) -> tuple[int, int, int]:
    """
    takes two strings a, b and the minimum overlap we're looking for
     returns where a starts, where b starts and the total overlap length
    Args:
        a (str):    a string
        b (str):    a string
        min_overlap (int):  the minimum overlap we're looking for
    Returns:
        tuple[int, int, int]:   where a starts, where b starts and the total overlap length
    """
    a_len, b_len = len(a), len(b)

    # we compare forwards and backwards at the same time, return the first largest match
    for b_len_to_min in range(b_len, min_overlap - 1, -1):
        # "1122|3344" "4455|" 4, "11223|344" "445|5" 3, "112233|44" "44|55" 2 MATCH

        if a[a_len - b_len_to_min:] == b[:b_len_to_min]:
            return a_len - b_len_to_min, 0, b_len_to_min

        # "0011" "1122|3344" 4, "0|111" "112|23344" 3, "00|11" "11|223344" 2 MATCH
        if b[b_len - b_len_to_min:] == a[:b_len_to_min]:
            return 0, b_len_to_min - 1, b_len_to_min

    # exit the for loop, we didn't find an overlap
    return 0, 0, 0


def assemble(seq1, seq2) -> str:
    """
    assemble two sequences
    :param seq1: a sequence
    :param seq2: a sequence
    :return: assembled string
    """

    overlap = best_overlap_with_min(seq1, seq2)
    if overlap[0] > overlap[1]:
        # seq1 starts before seq2
        return seq1 + seq2[overlap[2]:]
    else:
        # seq2 starts before seq1
        return seq2 + seq1[overlap[2]:]




def assemble_sequences_greedy(sequences, parallel= 1) -> str:
    """
    assemble all sequences into one
    Args: sequences (list[str]): a list of sequences
    Returns: str: the assembled sequence
    """

    #   pick firsts sequence, assemble it with the one with the smallest overlap, repeat
    while len(sequences) > 1:
        
        to_assemble = sequences[0]

        #IDX,
        best_overlap : (int, int) = (0, 0)

        #TODO optimize when parallel vs sequential (maybe sequential when under N sequences)
        if parallel and len(sequences) > 3000:
            best_overlap = parallel_find_best_pair(sequences)
        else:
            print("sequential, nb of seq: " + str(len(sequences)))
            sequences.pop(0)
            result = process_chunk_find_closest_pair(to_assemble, sequences, 0)
            best_overlap = (result[0], result[1])


        
        #we assemble the closest pair
        to_assemble = assemble(to_assemble, sequences[best_overlap[0]])
        sequences.pop(best_overlap[0])
        sequences.append(to_assemble)
    return sequences[0]



def process_chunk_find_closest_pair(sequence_to_assemble, sequence_chunk, chunk_number):
    #first is seq position, second  is the overlap length
    highest_overlap_with: (int, int) = (0,0)
    for i in range(len(sequence_chunk)):
        #    tuple[int, int, int]: where a starts, where b starts and the total overlap length
        current_tuple = best_overlap_with_min(sequence_to_assemble, sequence_chunk[i], highest_overlap_with[1])
        if current_tuple[2] > highest_overlap_with[1]:
            highest_overlap_with = (i, current_tuple[2])
    return (highest_overlap_with[0], highest_overlap_with[1], chunk_number)

    


"""
param: list of sequences (str)
return: best pair tuple (int, int)
"""
def parallel_find_best_pair(sequences, num_processes=None):
    if num_processes is None:
        #num_processes = multiprocessing.cpu_count()
        num_processes = 8

    print("parallel: nb of sequences: " + str(len(sequences)))
    sequence_to_assemble = sequences[0]
    sequences.pop(0)



    # Splitting the sequences into chunks
    chunk_size = len(sequences) // num_processes
    sequence_chunks = [sequences[i:i + chunk_size] for i in range(0, len(sequences), chunk_size)]



    # Prepare tuples of arguments for each chunk
    #chunk_data = [(sequence_to_assemble, chunk, chunk_id) for chunk_id, chunk in enumerate(sequence_chunks)]

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Map the sequence chunks to the worker processes
        results = pool.starmap(process_chunk_find_closest_pair, [(sequence_to_assemble, chunk, chunk_id) for chunk_id, chunk in enumerate(sequence_chunks)])

    highest_overlap_with = (0,0)
    for result in results:
        if result[1] > highest_overlap_with[1]:
            #result._3 is chunk number, we need to add (it*chunk size) to the index
            #result._1 is the index of the sequence in the chunk
            #result._2 is the overlap length
            highest_overlap_with = (result[0] + (result[2] * chunk_size), result[1])

    return highest_overlap_with






def time_assemble_sequences(path: str) -> str:
    append_string_to_file("\n___________\n")
    append_string_to_file(path)
    starttime = timeit.default_timer()
    result = assemble_sequences_greedy(fastaReader.fasta_to_list_str(path))
    append_string_to_file("\nThe time difference for custom is :" + str(timeit.default_timer() - starttime))
    append_string_to_file("\nlen of = " + str(len(result)))



def append_string_to_file(text_to_append, filename = "results/result3.txt"):
    # with open(filename, 'a') as file:
    #     file.write(text_to_append)
    print(text_to_append)

#
# print("greedy assembly on small reads fasta")
# ecoli_small_reads = fastaReader.fasta_to_list_str("fastas/Escherichia_coli_SMALL_READS.fasta")

# print(assemble_sequences_greedy(ecoli_small_reads))



if __name__ == "__main__":
    print("exo")
    print(assemble_sequences_greedy(["CTGTA", "ACCTG", "CCTGT"]))
    list_to_test_corona = ("fastas/corona/MN908947.3_fraction0000100_READS_MIXED.fasta",
                            "fastas/corona/MN908947.3_fraction0000500_READS_MIXED.fasta",
                            "fastas/corona/MN908947.3_fraction0001000_READS_MIXED.fasta",
                            "fastas/corona/MN908947.3_fraction0005000_READS_MIXED.fasta",
                            "fastas/corona/MN908947.3_fraction0010000_READS_MIXED.fasta",
                            "fastas/corona/MN908947.3_READS_MIXED.fasta")



    list_ecoli = ("fastas/ecoli/Escherichia_coli_fraction0000100_READS_MIXED.fasta", 
                "fastas/ecoli/Escherichia_coli_fraction0000500_READS_MIXED.fasta",
                "fastas/ecoli/Escherichia_coli_fraction0001000_READS_MIXED.fasta", 
                "fastas/ecoli/Escherichia_coli_fraction0005000_READS_MIXED.fasta",
                "fastas/ecoli/Escherichia_coli_fraction0010000_READS_MIXED.fasta",
                "fastas/ecoli/Escherichia_coli_fraction0050000_READS_MIXED.fasta",
                "fastas/ecoli/Escherichia_coli_fraction0100000_READS_MIXED.fasta",
                "fastas/ecoli/Escherichia_coli_fraction0500000_READS_MIXED.fasta")

    for fasta in list_ecoli:
        time_assemble_sequences(fasta)
