import re
from collections import OrderedDict
from urllib import request as req

# TODO: Improve function by testing different URL reader libraries against urllib.request()
def _grab_HTML_FASTA_data(protein_IDs):
    """ Grabs parsed FASTA data from HTML links and stores to dictionary. """
    pro_dict = dict()
    # Grabs and decodes protein data online and stores in dictionary
    for proID in protein_IDs:
        proURL = "http://www.uniprot.org/uniprot/{}.fasta".format(proID)
        url_data = req.urlopen(proURL)
        html = "".join(url_data.read().decode("utf-8").split("\n")[1:])
        pro_dict[proID] = html
    return pro_dict

def _format_output_motifs(motif_locations):
    """ Converts unstructured motif location data into structured output format. """
    output_format = str()
    # Structures motif positional data into final output format
    for motifID, motif_positions in motif_locations.items():
        output_format += "{}\n".format(motifID)
        output_format += " ".join([str(motif_pos) for motif_pos in motif_positions])
        output_format += "\n"
    return output_format

# TODO: Improve function by creating custom re.compile(), re.finditer(), OrderedDict() data structures
def get_motif_locations_across_FASTA(motifs_bank):
    """ Determines locations of motifs across FASTA values using algorithmic substring-searching. """
    regex, motif_locs = re.compile("(?=N[^P][ST][^P])"), OrderedDict()
    for seqID, sequence in motifs_bank.items():
        # Grabs overlapping RegEx matches across motif data
        matches = regex.finditer(sequence)
        # Generates motif positions data from RegEx generator
        if matches:
            motifs = list()
            for match in matches:
                motifs.append(match.start() + 1)
            #  Adds RegEx-matched motif positions in unstructured motifs dictionary
            if not motifs:
                continue
            else:
                motif_locs[seqID] = motifs
    return motif_locs

def main():
    # NOTE: Requires being in parent repo ('pwd' must return up to directory '/Rosalind_Bioinformatics/Bioinformatics_Stronghold')
    FILEPATHREAD = "rosalind_mprt.txt"
    FILEPATHWRITE = "rosalind_splc.txt"

    # Reads text data from raw dataset as single-line array of characters
    with open(FILEPATHREAD, "r") as fr:
        data = [value.strip() for value in fr.readlines()]

    # Sends FASTA dictionary to get_motifs() --> dict-value iterable must be WITHIN get_motifs()
    output_data = get_motif_locations_across_FASTA(_grab_HTML_FASTA_data(data))

    # Creates output file and writes appropriate response to file and notifies user
    with open(FILEPATHWRITE, "w") as fw:
        fw.write(_format_output_motifs(output_data)[:-1])

    return print("\nThe Protein Motifs dataset has been processed and the appropriate output has been saved to {}.\n".format(FILEPATHWRITE[2:]))

if __name__ == "__main__":
    main()